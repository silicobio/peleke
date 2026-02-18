
"""
Script to download a CIF structure from a given URL, convert it to PDB format, and score it using HADDOCK3.
"""

import argparse, os, json, re, subprocess
import logging


## Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

## Parse URL from command line arguments
parser = argparse.ArgumentParser(description="Scoring a PDB structure with HADDOCK3.")
parser.add_argument("--pdb_path", type=str, required=True, help="Path to the PDB file")
args = parser.parse_args()

pdb_path = args.pdb_path

def haddock3_score(pdb_path:str, output_dir: str) -> dict:
    """
    Runs the haddock3-score CLI tool on the given PDB file and extracts the scoring metrics.

    Parameters
    ----------
    pdb_path : str
        Path to the input PDB file.
    output_dir : str
        Directory to run the haddock3-score command in (where output files will be generated).

    Returns
    -------
    metrics : dict
        A dictionary containing the extracted scoring metrics.
    """
    try:
        ## Run haddock3-score CLI
        command = ["haddock3-score", "--full", pdb_path]
        sp_result = subprocess.run(command, cwd=output_dir, capture_output=True, text=True, check=True)

        ## Parse result
        metrics = {}

        ## Extract HADDOCK score
        match = re.search(r"HADDOCK-score \(emscoring\) = ([\-\d\.]+)", sp_result.stdout)
        if match:
            metrics["score"] = float(match.group(1))

        ## Extract individual energy terms
        matches = re.findall(r"(\w+)=([\-\d\.]+)", sp_result.stdout)
        for key, value in matches:
            metrics[key] = float(value)

        ## Calculate total score
        metrics["total"] = metrics["vdw"] + metrics["elec"]

        ## Remove air
        del metrics["air"]

        return metrics

    except subprocess.CalledProcessError as e:
        print("HADDOCK3 Error occurred:", e.stderr)
        return {}

  
if __name__ == "__main__":
    ## Make sure output directory exists
    output_dir = os.path.dirname(pdb_path)
    pdb_file_name = os.path.basename(pdb_path)
    os.makedirs(output_dir, exist_ok=True)

    ## Score with HADDOCK3
    logging.info(f"Scoring PDB structure with HADDOCK3: {pdb_path}")
    scores = haddock3_score(pdb_file_name, output_dir)
    
    ## Write scores dict to a json file
    if scores:
        output_file = f"{output_dir}/{pdb_file_name.split('.')[0]}_haddock3_scores.json"
        with open(output_file, "w") as f:
            json.dump(scores, f, indent=4)
        logging.info(f"\tHADDOCK3 scores saved to {output_file}")
    else:
        logging.warning(f"\tNo HADDOCK3 scores generated for {pdb_file_name}")