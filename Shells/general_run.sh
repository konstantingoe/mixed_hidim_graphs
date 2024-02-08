#!/bin/bash
#SBATCH -o ./outfiles/microbench.out
#SBATCH -D ./
#SBATCH -J microbenchmark
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=2000mb
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=72:00:00
#--------------------------------------

module load slurm_setup
module load r

cd /dss/dsshome1/lxc0A/ge43doh2/Projects/mixed_hidim_graphs

Rscript Functions/microbenchmark.R
