Abrar_Kurtosis.py: This script calculates mean, variance, kurtosis, sixth moment for different bin sizes but it can't tile the whole genome into the desired number of bins. Instead it tiles the AXT sequences into the desired number of bins.

Moment_plots.py: This script just plots the histograms from the outputs of the binned_moment_genome_ave.R script. 

export PYTHONPATH="/hpcshare/appsunit/MillerU/Python/3.10.1/lib/python3.10/site-packages/:$PYTHONPATH"
export PATH="/hpcshare/appsunit/MillerU/Python/3.10.1/bin:$PATH"
python3.10 "/bucket/.mabuya/MillerU/Abrar/Moment_val_bins/Abrar_Kurtosis.py"



#!/bin/bash
#SBATCH --job-name=python_Kurtosis
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=mdabrar.jahin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/
ml python/3.7.3
python