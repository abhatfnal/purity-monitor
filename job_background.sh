#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=ako.jamil@yale.edu
#SBATCH --partition=day
#SBATCH --job-name=cluster
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --time=10:00:00

python background.py
