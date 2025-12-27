#!/bin/bash
#SBATCH -J benchmark_run
#SBATCH -D ./MolSim-WS25-GroupD/
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
 
module load slurm_setup

./build/benchmarks/MolSimBench --benchmark_filter=Simulation/Full/
