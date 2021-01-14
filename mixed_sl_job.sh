#!/bin/bash
#SBATCH -o /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning/out
#SBATCH -D /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning
#SBATCH -J mixed_SL
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=72:00:00

module load slurm_setup

module load r

cd /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning

R –f server_run.R