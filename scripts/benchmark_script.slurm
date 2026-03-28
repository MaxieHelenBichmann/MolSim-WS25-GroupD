#!/bin/bash
#SBATCH -J benchmark_run
#SBATCH -D ./
#SBATCH -o ./results/benchmark/%j.out
#SBATCH -e ./results/benchmark/%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
 
module load slurm_setup

mkdir -p results/benchmark
./MolSim-WS25-GroupD/build/benchmarks/MolSimBench --benchmark_filter=Simulation/Full/
