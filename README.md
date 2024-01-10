# Covalent-Docking-Proximity-Analyzer
This Code is written by:
## Bikash Debnath

## Purpose:
Covalent docking is a very time and cost consuming process. For any covalent binding, first molecule needs to bind at the active site as intact form, then if warhead is close enough to reactive resudue, then it forms covalent bond. On this principle, I used normal docking as docking, then used those docking poses in this analysis to find out warheads of which compounds (within top 3 poses) are acually docking proximity to reactive residue. Then filter those compounds and redock using covalent docking using choice of covalent docking methods.

This code has been optimized in GCP slurm.
## Requirements 
1. RDKit
2. OpenBabel
3. Prolif
4. MDAnalysis

## create a separate virtual environment
conda create -n prolif
## activate it
conda activate prolif
conda install -c conda-forge prolif

## Create multiple folders: poses, input_txt, log, subtract, output
## prepare protein (protein.pdb) and keep it in the working directory, you can use a molecule as reference (ref.sdf) in the working directory.

## Keep all individual docked pdbqt files from vina docking in the poses folder; You can use other formats also

## To parallel processing split folder into multiple folders

$python folder_splitter1.py poses 50

## Follow these two steps to convert poses into sdf file with all Hydrogens

$for d in poses/*/; do (cd "$d" && obabel *.pdbqt -O ligands-dH.mol2 -d ---errorlevel 0); done;

$for d in poses/*/; do (cd "$d" && obabel ligands-dH.mol2 -O ligands-hH.sdf -p 7.4 ---errorlevel 0); done;
 
## prepare input file to run bath array jobs
$python write_input_txt.py poses

## subtract covalent linker from sdf file. Basically, Here we want only covalent linker from the docked pose to calculte proximity between covalent linker and Cystein residue of the protein. Based on demand you can add warheads in "warheads_subtract.py": core and sub_core smarts. I tried to cover few common warheads. 
sbatch --array=1-50 cov_subtract.sh
## Activate prolif environment, if not
$conda activate prolif
## Calculate proximity
### Distance is set at 4 Angstrom, You can modify by going the "prolif_cov.py" at "class CustomHydrophobic(plf.interactions.Hydrophobic)" Generally, you can set the distance to any value between 3.5 - 5 angstrom. Change Residue number based on your protein in this file.
$sbatch --array=1-50 slurmscriptfinal.sh
## Filter the poses which are close proximity of 4 Angstrom
$sbatch --array=1-50 pose_filter.sh
