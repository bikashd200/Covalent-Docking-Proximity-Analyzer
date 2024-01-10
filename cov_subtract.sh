#!/usr/bin/env bash

# Slurm Settings
###############################################################################

#SBATCH --job-name=Prolif
#SBATCH --mail-user=bdebnath@championsoncology.com
#SBATCH --mail-type=fail
#SBATCH --time=03-12:00:00
#SBATCH --nodes=1
#SBATCH --partition=p2
#SBATCH --output=./log/subtract_%A_%a.out           # File to which standard out will be written
#SBATCH --error=./log/subtract_%A_%a.err           # File to which standard  err will be written
#SBATCH --signal=10@300


# Do something else.
echo Started at
#python times.py
date

python warheads_subtract.py input-$SLURM_ARRAY_TASK_ID.txt

echo Completed at
date
#python times.py

