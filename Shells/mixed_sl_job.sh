#!/bin/bash
#SBATCH -o /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning/mixed_binary_tiny.out
#SBATCH -D /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning
#SBATCH -J "mixed_binary"
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=konstantin.goebler@tum.de
#SBATCH --export=NONE
#SBATCH --time=200:00:00

module load slurm_setup

R_LIBS_USER="/dss/dsshome1/lxc0A/ge43doh2/R/x86_64-pc-linux-gnu-library/3.6" 
export R_LIBS_USER
R_LIBS_SITE="/dss/dsshome1/lrz/sys/spack/release/19.2/opt/x86_avx2/r/3.6.0-gcc-eujetj3/rlib/R/library"
export R_LIBS_SITE

module unload intel-mpi
module unload intel
module load gcc
module load r

cd /dss/dsshome1/lxc0A/ge43doh2/mixed/mixed_structure_learning

Rscript binary_run.R
