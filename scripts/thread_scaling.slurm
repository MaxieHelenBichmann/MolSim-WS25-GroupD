#!/bin/bash
#SBATCH -J thread_scaling
#SBATCH -D ./
#SBATCH -o ./results/benchmark/%j.out
#SBATCH -e ./results/benchmark/%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=01:00:00

module load slurm_setup

mkdir -p results/benchmark
./MolSim-WS25-GroupD/build/benchmarks/MolSimBench --benchmark_filter=ThreadScaling --benchmark_out=results/benchmark/thread_scaling_${SLURM_JOB_ID}.json --benchmark_out_format=json
