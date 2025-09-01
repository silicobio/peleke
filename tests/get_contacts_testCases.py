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

def get_epitope_residues_pandaprot(pdb_file, chain_set1, chain_set2):
    """
    Extracts epitope residues from a PDB file using PandaProt.
    Args:
        pdb_file: str, path to the PDB file (can be gzipped).
        chain_set1: list of str, first set of chain identifiers (e.g., ['H', 'L']).
        chain_set2: list of str, second set of chain identifiers (e.g., ['A', 'B', 'C']).
    Returns:
        tuple: (epitope_residues_set1, epitope_residues_set2) where each is a list of str in the format 'A:ARG 176',
    """
    try:
        logger.info("1. Running PandaProt analysis...")
        chains = chain_set1 + chain_set2
        ## Make the PandaProt stuff quiet
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                analyzer = PandaProt(pdb_file, chains=chains)
                interactions = analyzer.map_interactions()


        epitope_residues_set1 = []
        epitope_residues_set2 = []
        relevant_interactions = []
        for interaction_type, interactions_list in interactions.items():
            for interaction in interactions_list:
                chain1 = interaction.get('chain1', interaction.get('donor_chain', ''))
                chain2 = interaction.get('chain2', interaction.get('acceptor_chain', ''))
                res1 = interaction.get('residue1', interaction.get('donor_residue', ''))
                res2 = interaction.get('residue2', interaction.get('acceptor_residue', ''))
                ## Consider interactions between the two chain sets
                if (chain1 in chain_set1 and chain2 in chain_set2):
                    epitope_residues_set1.append(f"{chain1}:{res1}")
                    epitope_residues_set2.append(f"{chain2}:{res2}")
                    relevant_interactions.append((interaction_type, interaction))
                elif (chain1 in chain_set2 and chain2 in chain_set1):
                    epitope_residues_set1.append(f"{chain2}:{res2}")
                    epitope_residues_set2.append(f"{chain1}:{res1}")
                    relevant_interactions.append((interaction_type, interaction))
                ## Ignores intra-set interactions

        return sorted(set(epitope_residues_set1)), sorted(set(epitope_residues_set2))
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return [], []

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


def find_contacts(pdb_id, pdb_file, chain_set1, chain_set2, target_chains, target_seqs, output_file):
    """
    Extracts epitope residues from a PDB file using PandaProt and highlights them in the target chain sequences.
    Args:
        pdb_id: str, PDB identifier.
        pdb_file: str, path to the PDB file (can be gzipped).
        chain_set1: list of str, first set of chain identifiers (e.g., ['H', 'L']).
        chain_set2: list of str, second set of chain identifiers (e.g., ['A', 'B', 'C']).
        target_chains: list of str, identifiers for chains to analyze epitopes (usually chain_set2).
        target_seqs: list of str, sequences corresponding to target_chains.
        output_file: str, path to output CSV file.
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
    required_chains = set(chain_set1 + chain_set2)
    
    ## Only use chains that are present
    chain_set1 = [c for c in chain_set1 if c in available_chains]
    chain_set2 = [c for c in chain_set2 if c in available_chains]
    target_chains = [c for c in target_chains if c in available_chains]

    if not chain_set1 or not chain_set2 or not target_chains:
        return logger.warning(f"Skipping {pdb_file}: Required chains ({required_chains}) not found. Available: ({available_chains})")
    
    # pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    residues_set1, residues_set2 = get_epitope_residues_pandaprot(pdb_file, chain_set1, chain_set2)

    chain_list, seq_list, res_list = [], [], []
    
    # Determine which residue set to use based on target chains
    # If target_chains are in chain_set2, use residues_set2; otherwise use residues_set1
    if all(chain in chain_set2 for chain in target_chains):
        target_residues = residues_set2
    else:
        target_residues = residues_set1
    
    for i, target_chain in enumerate(target_chains):
        # If target_seqs is None, extract sequence from PDB
        if target_seqs is None:
            atom_df = pdb_df.df['ATOM']
            chain_df = atom_df[atom_df['chain_id'] == target_chain]
            residues = chain_df[['residue_number', 'residue_name']].drop_duplicates()
            aa_map = {
                'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
                'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
            }
            target_sequence = ''.join([aa_map.get(resn, 'X') for _, resn in residues.values])
        else:
            target_sequence = target_seqs[i] if i < len(target_seqs) else None
        if target_sequence and target_sequence != 'nan':
            try:
                resnum_to_idx = build_resnum_to_seq_idx_map(pdb_df, target_chain)
                highlighted_seq = highlight_epitope_in_sequence(target_sequence, target_chain, target_residues, resnum_to_idx)
            except Exception as e:
                print(f"Error mapping for {pdb_id} chain {target_chain}: {e}")
                highlighted_seq = None
        else:
            highlighted_seq = None
        chain_list.append(target_chain)
        seq_list.append(highlighted_seq if highlighted_seq else "")
        ## Only include residues for this chain
        chain_residues = [r for r in target_residues if r.startswith(f"{target_chain}:")]
        res_list.append('|'.join(chain_residues))

        combined_results = [{
            'pdb_id': pdb_id,
            'target_chain_ids': '|'.join(chain_list),
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
    parser.add_argument("--chain_set1", type=str, required=True, help="List of |-delimited chain identifiers for first set (e.g., 'H|L').")
    parser.add_argument("--chain_set2", type=str, required=True, help="List of |-delimited chain identifiers for second set (e.g., 'A|B|C').")
    parser.add_argument("--target_chains", type=str, required=True, help="List of |-delimited target chain identifiers for epitope analysis (e.g., 'A|B|C').")
    # No longer need --target_seqs; sequences will be extracted from PDB
    parser.add_argument("--output_file", type=str, required=True, help="Output .csv file to save results.")
    
    # Add backward compatibility arguments (optional)
    parser.add_argument("--h_chain_id", type=str, help="Heavy chain identifier (for backward compatibility).")
    parser.add_argument("--l_chain_id", type=str, help="Light chain identifier (for backward compatibility).")
    parser.add_argument("--antigen_ids", type=str, help="Antigen chain identifiers (for backward compatibility).")
    parser.add_argument("--antigen_seqs", type=str, help="Antigen sequences (for backward compatibility).")
    
    args = parser.parse_args()

    # Handle backward compatibility
    if args.h_chain_id and args.l_chain_id and args.antigen_ids and args.antigen_seqs:
        logger.info("Using backward compatibility mode (antibody-antigen analysis)")
        chain_set1 = [args.h_chain_id, args.l_chain_id]
        chain_set2 = args.antigen_ids.split('|')
        target_chains = chain_set2
        target_seqs = args.antigen_seqs.split('|')
    else:
        chain_set1 = args.chain_set1.split('|')
        chain_set2 = args.chain_set2.split('|')
        target_chains = args.target_chains.split('|')
        target_seqs = None  # Will be extracted from PDB

    logger.info(f"Processing {args.pdb_file} for chain set 1: {chain_set1}, chain set 2: {chain_set2}, target chains: {target_chains}")

    pdb_file = args.pdb_file
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    output_file = args.output_file

    ## Check if output file exists. If so, skip processing
    if os.path.exists(output_file):
        logger.warning(f"Output file {output_file} already exists. Skipping processing for {pdb_id}.")
    else:
        logger.info(f"Starting contact extraction for {pdb_id}...")
        find_contacts(pdb_id, pdb_file, chain_set1, chain_set2, target_chains, target_seqs, output_file)