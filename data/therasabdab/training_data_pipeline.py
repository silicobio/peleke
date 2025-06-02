import os
import subprocess
import requests
import pandas as pd
from Bio.SeqUtils import seq1

# --- Configuration ---
source_file = "TheraSAbDab_SeqStruc_OnlineDownload.csv"
pdb_dir = "/mnt/data/pdb_files"
os.makedirs(pdb_dir, exist_ok=True)

# --- Load Dataset ---
df = pd.read_csv(source_file)
structure_cols = ['100% SI Structure', '99% SI Structure', '95-98% SI Structure']

def get_pdb_id(row):
    for col in structure_cols:
        val = row[col]
        if pd.notna(val) and val.lower() != 'na':
            return val.split(';')[0].strip().lower()
    return None

df['PDB_ID'] = df.apply(get_pdb_id, axis=1)
df = df[df['PDB_ID'].notna()]

# --- Helper: Download PDB ---
def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(out_path):
        r = requests.get(url)
        r.raise_for_status()
        with open(out_path, 'w') as f:
            f.write(r.text)
    return out_path

# --- Helper: Run PyMOL to Find Polar Contacts ---
def run_pymol_find_contacts(input_path, output_path):
    script_path = os.path.join(pdb_dir, "find_polar_contacts_script.py")
    script_content = f"""
from pymol import cmd

cmd.load('{input_path}', 'complex')
cmd.remove('hydro')
cmd.remove('not (chain A+B)')

for chain in ['A', 'B']:
    model = cmd.get_model(f'chain {{chain}}')
    resvs = [atom.resv for atom in model.atom]
    if not resvs:
        continue
    min_resv = min(resvs)
    cmd.alter(f'chain {{chain}}', f'resv = resv - {{min_resv}} + 1')

cmd.rebuild()
cmd.save('{input_path}', 'complex')

cmd.select('ab', 'chain A')
cmd.select('ag', 'chain B')

pairs = cmd.find_pairs('ab and name N+O', 'ag and name N+O', cutoff=3.5)

contact_residues = set()
for ab_atom, ag_atom in pairs:
    contact_residues.add(ag_atom[3])

with open('{output_path}', 'w') as f:
    f.write(','.join(str(r) for r in sorted(contact_residues)))

cmd.quit()
"""
    with open(script_path, 'w') as f:
        f.write(script_content)
    subprocess.run(["pymol", "-cq", script_path], check=True)

# --- Extract Sequence from PDB File with Brackets ---
def extract_sequence_with_brackets(pdb_path, chain, contact_residues):
    seq = ""
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain:
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    continue
                res_name = line[17:20].strip()
                key = (res_num, res_name)
                if key in seen:
                    continue
                seen.add(key)
                try:
                    aa = seq1(res_name)
                except KeyError:
                    aa = 'X'  # unknown residue
                if res_num in contact_residues:
                    seq += f"[{aa}]"
                else:
                    seq += aa
    return seq

# --- Process Entries ---
output_rows = []
for _, row in df.iterrows():
    pdb_id_full = row['PDB_ID']
    pdb_code = pdb_id_full.split(':')[0]
    chain_ids = pdb_id_full.split(':')[1] if ':' in pdb_id_full else ''

    try:
        pdb_path = download_pdb(pdb_code)
        mod_pdb_path = os.path.join(pdb_dir, f"mod_{pdb_code}.pdb")
        contacts_out_path = os.path.join(pdb_dir, f"contacts_{pdb_code}.txt")

        # Create modified PDB with antibody chains -> A, antigen chains -> B
        with open(pdb_path) as fin, open(mod_pdb_path, 'w') as fout:
            for line in fin:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    if chain_id in chain_ids:
                        line = line[:21] + 'A' + line[22:]
                    else:
                        line = line[:21] + 'B' + line[22:]
                fout.write(line)

        run_pymol_find_contacts(mod_pdb_path, contacts_out_path)

        with open(contacts_out_path) as f:
            contacts_str = f.read().strip()
            contacts = list(map(int, contacts_str.split(','))) if contacts_str else []

        antigen_seq = extract_sequence_with_brackets(mod_pdb_path, 'B', set(contacts))

        heavy = row['HeavySequence(ifbispecific)'] if pd.notna(row['HeavySequence(ifbispecific)']) else row['HeavySequence']
        light = row['LightSequence(ifbispecific)'] if pd.notna(row['LightSequence(ifbispecific)']) else row['LightSequence']
        if pd.isna(heavy) or pd.isna(light):
            continue
        antibody_seq = f"{heavy}|{light}"

        output_rows.append({"antigen": antigen_seq, "antibody": antibody_seq})

    except Exception as e:
        print(f"Failed on {pdb_code}: {e}")
        continue

# --- Save Output ---
output_df = pd.DataFrame(output_rows)
output_path = "/mnt/data/antibody_training_dataset.csv"
output_df.to_csv(output_path, index=False)
print(f"Saved training data to {output_path}")
