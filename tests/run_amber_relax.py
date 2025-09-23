import os
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer
import logging
import torch

"""
Minimize (relax the side chains of) all PDB files in a given directory using OpenMM.
Inputs:
    - input_dir:str The input directory of PDB files to be minimized
    - output_dir:str The output directory where minimized PDB files will be saved
Outputs:
    - Minimized PDB files saved in the output directory with the same names as the input files.
"""

## Set up logging
logging.basicConfig(level=logging.INFO)

def minimize_pdb(input_pdb, output_pdb):
    ## Load and prepare structure
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    ## Create system using AMBER force field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )

    ## Integrator (not used for MD, just needed)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    if torch.cuda.is_available():
        logging.info("Using GPU for OpenMM.")
        platform = Platform.getPlatformByName('CUDA')
        properties = {"CudaDeviceIndex": "0", "CudaPrecision": "mixed"}
        # platform = Platform.getPlatformByName('OpenCL')
        # properties = {"OpenCLPrecision": "mixed"}
    else:
        logging.info("Using CPU for OpenMM.")
        platform = Platform.getPlatformByName('CPU')
        properties = {}
    simulation = Simulation(fixer.topology, system, integrator, platform, properties)
    simulation.context.setPositions(fixer.positions)

    ## Energy minimization
    simulation.minimizeEnergy()
    positions = simulation.context.getState(getPositions=True).getPositions()

    ## Save output
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, positions, f)

def main():
    parser = argparse.ArgumentParser(description="Minimize PDBs using OpenMM.")
    parser.add_argument("input_dir", help="Directory with input PDB files")
    parser.add_argument("output_dir", help="Directory to save minimized PDB files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")

    pdb_files = [f for f in os.listdir(args.input_dir) if f.endswith(".pdb")]
    if not pdb_files:
        logging.warning("No PDB files found in input directory.")
        return

    for pdb_file in pdb_files:
        input_path = os.path.join(args.input_dir, pdb_file)
        output_path = os.path.join(args.output_dir, pdb_file)
        logging.info(f"Minimizing: {pdb_file}")
        try:
            minimize_pdb(input_path, output_path)
        except Exception as e:
            logging.error(f"Failed to minimize {pdb_file}: {e}")

if __name__ == "__main__":
    main()
