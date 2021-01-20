#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ako.jamil@yale.edu
#SBATCH --partition=day
#SBATCH --job-name=YLXPM
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=10:00:00


python Background.py
