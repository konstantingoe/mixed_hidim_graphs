#!/bin/bash
#SBATCH -o ./outfiles/general.out
#SBATCH -D ./
#SBATCH -J general_run
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=2600000mb
#SBATCH --cpus-per-task=50
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=18:00:00
#--------------------------------------

module load slurm_setup
module load r

cd /dss/dsshome1/lxc0A/ge43doh2/Projects/mixed_hidim_graphs

Rscript /Functions/new_sims.R
