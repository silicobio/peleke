import re, os
import subprocess
import pandas as pd

## Find all .pdb files in a directory
structures_dir = '/mnt/tests/structures/predicted_complexes'
pdb_files = [f for f in os.listdir(structures_dir) if f.endswith('.pdb')]

def haddock3_score(pdb_path:str) -> dict:
  try:
    ## Run haddock3-score CLI
    command = ["haddock3-score", "--full", pdb_path]
    sp_result = subprocess.run(command, capture_output=True, text=True, check=True)

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
  

## Pandas dataframe to store results
affinity_df = pd.DataFrame(columns=['seq_id', 'score', 'total', 'vdw', 'elec', 'desolv', 'bsa'])
  
## Run haddock3_score on each pdb file and append results to dataframe
for pdb_file in pdb_files:
    pdb_path = os.path.join(structures_dir, pdb_file)
    print(f"Scoring {pdb_file}")
    scores = haddock3_score(pdb_path)
    if scores:
        scores['seq_id'] = pdb_file.replace('.pdb','')
        affinity_df_row = pd.DataFrame([scores])
        affinity_df = pd.concat([affinity_df, affinity_df_row], ignore_index=True)


## Save results to CSV
output_file = 'binding_affinity_results.csv'
affinity_df.to_csv(output_file, index=False)