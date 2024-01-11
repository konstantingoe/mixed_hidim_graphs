#!/bin/bash
#SBATCH -o $HOME/Projects/mixed_hidim_graphs/outfiles/serial.out
#SBATCH -D $HOME/Projects/mixed_hidim_graphs
#SBATCH -J "mixed_serial"
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_long
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=8-00:00:00

module load slurm_setup
module load r

cd $HOME/Projects/mixed_hidim_graphs/Functions

Rscript new_sims.R
