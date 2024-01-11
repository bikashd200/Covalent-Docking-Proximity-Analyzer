#!/usr/bin/env bash

# Slurm Settings
###############################################################################

#SBATCH --job-name=Prolif
#SBATCH --mail-user=user@mail.com
#SBATCH --mail-type=fail
#SBATCH --time=03-12:00:00
#SBATCH --nodes=1
#SBATCH --partition=p1
#SBATCH --output=./log/prolif_%A_%a.out           # File to which standard out will be written
#SBATCH --error=./log/prolif_%A_%a.err           # File to which standard  err will be written
#SBATCH --signal=10@300


# Do something else.
echo Started at
#python times.py
date

python prolif_cov.py input-$SLURM_ARRAY_TASK_ID.txt

echo Completed at
date
#python times.py

