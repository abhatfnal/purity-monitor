#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ako.jamil@yale.edu
#SBATCH --partition=day
#SBATCH --job-name=plotdataall
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --time=1:00:00
loadconda 
python PlotAllData.py 
