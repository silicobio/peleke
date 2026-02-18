import os, glob

pdb_files = glob.glob("../**/*.pdb", recursive=True)

for pdb_file in pdb_files:
    ## Call the score script for each PDB file
    command = f"sbatch submit.sh {pdb_file}"
    id = os.path.basename(pdb_file).split('.')[0]
    print(f"Submitting HADDOCK3 scoring job for {id}")
    os.system(command)