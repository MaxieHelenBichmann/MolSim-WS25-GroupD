#!/bin/bash
#SBATCH -J intel-profile
#SBATCH -D ./
#SBATCH -o ./results/intel/%j.out
#SBATCH -e ./results/intel/%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
 
module load slurm_setup
module load intel-oneapi-vtune/2024.2.1 
module load intel-oneapi-advisor/2024.2.1

mkdir -p results/intel
make -C MolSim-WS25-GroupD/build intel-profile-all
