from pandaprot import PandaProt
from biopandas.pdb import PandasPdb
import os, re
import pandas as pd
import tempfile
import shutil
import gzip
import argparse
import logging
from contextlib import redirect_stdout

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_epitope_residues_pandaprot(pdb_file, h_chain_id, l_chain_id, antigen_ids):
    """
    Extracts epitope residues from a PDB file using PandaProt.
    Args:
        pdb_file: str, path to the PDB file (can be gzipped).
        h_chain_id: str, heavy chain identifier (e.g., 'H').
        l_chain_id: str, light chain identifier (e.g., 'L').
        antigen_ids: list of str, identifiers for antigen chains (e.g., ['A', 'B', 'C']).
    Returns:
        list of str: epitope residues in the format 'A:ARG 176',
    """
    try:
        logger.info("1. Running PandaProt analysis...")
        chains = [h_chain_id, l_chain_id] + antigen_ids
        ## Make the PandaProt stuff quiet
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                analyzer = PandaProt(pdb_file, chains=chains)
                interactions = analyzer.map_interactions()


        epitope_residues = []
        relevant_interactions = []
        for interaction_type, interactions_list in interactions.items():
            for interaction in interactions_list:
                chain1 = interaction.get('chain1', interaction.get('donor_chain', ''))
                chain2 = interaction.get('chain2', interaction.get('acceptor_chain', ''))
                res1 = interaction.get('residue1', interaction.get('donor_residue', ''))
                res2 = interaction.get('residue2', interaction.get('acceptor_residue', ''))
                ## Only consider antigen-antibody interactions
                for antigen_id in antigen_ids:
                    if (
                        (chain1 == antigen_id and chain2 in [h_chain_id, l_chain_id]) or
                        (chain2 == antigen_id and chain1 in [h_chain_id, l_chain_id])
                    ):
                        epitope_residues.append(f"{antigen_id}:{res1 if chain1 == antigen_id else res2}")
                        relevant_interactions.append((interaction_type, interaction))
                    ## Ignores antigen-antigen interactions

        return sorted(set(epitope_residues))
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return []

def build_resnum_to_seq_idx_map(pdb_df: PandasPdb, chain_id: str) -> dict:
    """
    Returns a dict mapping PDB residue numbers (as int) to sequence indices (1-based) for a given chain.
    Args:
        pdb_df: PandasPdb, BioPandas DataFrame of the PDB file.
        chain_id: str, chain identifier (e.g., 'A')
    Returns:
        dict: mapping of PDB residue numbers to sequence indices.
    """
    logger.info("2. Running BioPandas renumbering...")
    atom_df = pdb_df.df['ATOM']
    ## Only keep rows for the specified chain
    chain_df = atom_df[atom_df['chain_id'] == chain_id]
    ## Get unique residue numbers in order of appearance
    residues = chain_df[['residue_number', 'residue_name']].drop_duplicates()
    resnum_to_idx = {}
    for idx, (resnum, _) in enumerate(residues.values, 1):  # 1-based index
        resnum_to_idx[int(resnum)] = idx
    return resnum_to_idx

def highlight_epitope_in_sequence(sequence, chain_id, epitope_residues, resnum_to_idx):
    """
    Places brackets around epitope residues in the antigen sequence.
    Args:
        sequence: str, full antigen sequence (1-letter code)
        chain_id: str, chain identifier (e.g., 'A')
        epitope_residues: list of str, e.g., ['A:ARG 176', ...]
        resnum_to_idx: dict mapping PDB residue numbers to sequence indices (1-based)
    Returns:
        str: sequence with epitope residues highlighted in square brackets.
    """
    logger.info("3. Highlighting epitope residues in sequence...")
    import re
    pattern = re.compile(rf"^{chain_id}:(?:\w+)\s*(\d+)$")
    epitope_seq_indices = set()
    for res in epitope_residues:
        m = pattern.match(res)
        if m:
            pdb_resnum = int(m.group(1))
            seq_idx = resnum_to_idx.get(pdb_resnum)
            if seq_idx:
                epitope_seq_indices.add(seq_idx)
    highlighted_seq = ""
    for i, aa in enumerate(sequence, 1):
        if i in epitope_seq_indices:
            highlighted_seq += f"[{aa}]"
        else:
            highlighted_seq += aa
    return highlighted_seq


def find_contacts(pdb_id, pdb_file, h_chain_id, l_chain_id, antigen_ids, antigen_seqs, output_file):
    """
    Extracts epitope residues from a PDB file using PandaProt and highlights them in the antigen sequence.
    Args:
        pdb_file: str, path to the PDB file (can be gzipped).
        h_chain_id: str, heavy chain identifier (e.g., 'H').
        l_chain_id: str, light chain identifier (e.g., 'L').
        antigen_ids: list of str, identifiers for antigen chains (e.g., ['A', 'B', 'C']).
    Returns:
        str: message indicating processing status and results.
    """
    ## If .gz, unzip to a temp file
    if pdb_file.endswith('.gz'):
        with tempfile.NamedTemporaryFile(suffix='.pdb', mode='wb', delete=False) as temp_file:
        ## Open the .gz file in binary read mode
            with gzip.open(pdb_file, 'rb') as f:
                # Copy the unzipped content from the .gz file to the temporary file
                shutil.copyfileobj(f, temp_file)

        pdb_file = temp_file.name
        logger.info(f"Unzipped input PDB to temporary file {temp_file.name}")
    else:
        pdb_file = pdb_file

    pdb_df = PandasPdb().read_pdb(pdb_file)
    
    ## Get available chains and check if required chains are present
    available_chains = set(str(c).strip() for c in pdb_df.df['ATOM']['chain_id'].unique())
    required_chains = {h_chain_id, l_chain_id} | set(antigen_ids)
    
    ## Only use chains that are present
    h_chain_id = h_chain_id if h_chain_id in available_chains else None
    l_chain_id = l_chain_id if l_chain_id in available_chains else None
    antigen_ids = [c for c in antigen_ids if c in available_chains]

    if not h_chain_id or not l_chain_id or not antigen_ids:
        return logger.warning(f"Skipping {pdb_file}: Required chains ({required_chains}) not found. Available: ({available_chains})")
    
    # pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    residues = get_epitope_residues_pandaprot(pdb_file, h_chain_id, l_chain_id, antigen_ids)

    chain_list, seq_list, res_list = [], [], []
    for i, antigen_chain in enumerate(antigen_ids):
        antigen_sequence = antigen_seqs[i] if i < len(antigen_seqs) else None
        if antigen_sequence and antigen_sequence != 'nan':
            try:
                resnum_to_idx = build_resnum_to_seq_idx_map(pdb_df, antigen_chain)
                highlighted_seq = highlight_epitope_in_sequence(antigen_sequence, antigen_chain, residues, resnum_to_idx)
            except Exception as e:
                print(f"Error mapping for {pdb_id} chain {antigen_chain}: {e}")
                highlighted_seq = None
        else:
            highlighted_seq = None
        chain_list.append(antigen_chain)
        seq_list.append(highlighted_seq if highlighted_seq else "")
        ## Only include residues for this chain
        chain_residues = [r for r in residues if r.startswith(f"{antigen_chain}:")]
        res_list.append('|'.join(chain_residues))
        # combined_results.append()

        combined_results = [{
            'pdb_id': pdb_id,
            'antigen_ids': '|'.join(chain_list),
            'highlighted_epitope_seqs': '|'.join(seq_list),
            'epitope_residues': '|'.join(res_list)
        }]

    ## Create DataFrame and merge with original
    highlight_df = pd.DataFrame(combined_results)

    ## Write to CSV if output file is specified
    highlight_df.to_csv(output_file, index=False)
    return logger.info(f"DONE: {pdb_id} processed. Results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract epitope residues from PDB files using PandaProt.")
    parser.add_argument("--pdb_file", type=str, required=True, help="Path to the PDB file (can be gzipped).")
    parser.add_argument("--h_chain_id", type=str, required=True, help="Heavy chain identifier (e.g., 'H').")
    parser.add_argument("--l_chain_id", type=str, required=True, help="Light chain identifier (e.g., 'L').")
    parser.add_argument("--antigen_ids", type=str, required=True, help="List of |-delimited antigen chain identifiers (e.g., 'A|B|C').")
    parser.add_argument("--antigen_seqs", type=str, required=True, help="List of |-delimited antigen sequences.")
    parser.add_argument("--output_file", type=str, required=True, help="Output .csv file to save results.")
    
    args = parser.parse_args()

    logger.info(f"Processing {args.pdb_file} for antibody chains {args.h_chain_id}, {args.l_chain_id} with antigen chain(s) {args.antigen_ids}")

    pdb_file = args.pdb_file
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    h_chain_id = args.h_chain_id
    l_chain_id = args.l_chain_id
    antigen_ids = args.antigen_ids.split('|')
    antigen_seqs = args.antigen_seqs.split('|')
    output_file = args.output_file

    ## Check if output file exists. If so, skip processing
    if os.path.exists(output_file):
        logger.warning(f"Output file {output_file} already exists. Skipping processing for {pdb_id}.")
    else:
        logger.info(f"Starting contact extraction for {pdb_id}...")
        find_contacts(pdb_id, pdb_file, h_chain_id, l_chain_id, antigen_ids, antigen_seqs, output_file)