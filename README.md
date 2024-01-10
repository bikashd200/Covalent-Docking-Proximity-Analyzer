# Covalent-Docking-Proximity-Analyzer
This Code is written by:
## Bikash Debnath
This code has been optimized in GCP slurm
## Requirements 
1. RDKit
2. OpenBabel
3. Prolif

# create a separate virtual environment
conda create -n prolif
# activate it
conda activate prolif
conda install -c conda-forge prolif

# Keep all individual docked pdbqt files from vina docking in the poses folder; You can use other formats also

# To parallel processing split folder into multiple folders

$python folder_splitter1.py poses 50

## Follow these two steps to convert poses into sdf file with all Hydrogens

$for d in poses/*/; do (cd "$d" && obabel *.pdbqt -O ligands-dH.mol2 -d ---errorlevel 0); done;

$for d in poses/*/; do (cd "$d" && obabel ligands-dH.mol2 -O ligands-hH.sdf -p 7.4 ---errorlevel 0); done;
 
## prepare input file to run bath array jobs
$python write_input_txt.py poses

## subtract covalent linker from sdf file. Basically, Here we want only covalent linker from the docked pose to calculte proximity between covalent linker and Cystein residue of the protein
sbatch --array=1-50 cov_subtract.sh
## Activate prolif environment, if not
$conda activate prolif
## Calculate proximity
$sbatch --array=1-50 slurmscriptfinal.sh
## Filter the poses which are close proximity of 4 Angstrom
$sbatch --array=1-50 pose_filter.sh
