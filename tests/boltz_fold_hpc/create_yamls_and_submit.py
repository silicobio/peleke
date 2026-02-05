import pandas as pd
import os
import yaml

## Read experiment data from Excel file
experiment_df = pd.read_excel("../test_cases.xlsx", sheet_name="experiments")

## Define function to create Boltz-2 YAML configuration for multimer folding
def create_boltz_yaml(h_chain: str, l_chain: str, a_chain: str, output_path: str):
    config_json = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "H",
                    "sequence": h_chain,
                }
            },
            {
                "protein": {
                    "id": "L",
                    "sequence": l_chain,
                }
            },
            {
                "protein": {
                    "id": "A",
                    "sequence": a_chain,
                }
            }
        ]
        }
    
    ## Make directories if they don't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    ## Write YAML file
    print(f"Creating YAML for {seq_id} at {yaml_path}")
    with open(output_path, 'w') as f:
        yaml.dump(config_json, f)


## Define a function for running Boltz-2 folding from the command line
def run_boltz_folding(yaml_path: str, accelerator: str = 'gpu'):
    # boltz_command = f"boltz predict {yaml_path} --write_full_pae --accelerator {accelerator}"
    boltz_command = f"sbatch submit.sh {yaml_path}"
    print(f"Running command: {boltz_command}")
    os.system(boltz_command)


for index, row in experiment_df.iterrows():
    seq_id = row['seq_id']
    h_chain = row['h_chain']
    l_chain = row['l_chain']
    a_chain = row['a_chain']
    yaml_path = f"./{seq_id}.yaml"

    create_boltz_yaml(h_chain, l_chain, a_chain, yaml_path)
    run_boltz_folding(yaml_path, accelerator='gpu')