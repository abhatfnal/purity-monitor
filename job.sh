#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=ako.jamil@yale.edu
#SBATCH --partition=day
#SBATCH --job-name=cluster
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=0:01:00

python read.py -f /scratch/hep/david_moore/zl423/PurityData/20180606/1d4grids0.1in/1bar/field_80Vcm_ratio1.55/
