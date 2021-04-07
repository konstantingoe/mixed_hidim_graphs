#!/bin/bash
#SBATCH -o /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning/level_sensitivity.out
#SBATCH -D /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning
#SBATCH -J "level_sens"
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --cpus-per-task=18
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=12:00:00

module load slurm_setup

R_LIBS_USER="/dss/dsshome1/lxc0A/ge43doh2/R/debian10/3.6/"
export R_LIBS_USER

module load r

cd /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning

Rscript level_sensitivity.R
