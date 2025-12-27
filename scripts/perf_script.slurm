#!/bin/bash
#SBATCH -J perf_run
#SBATCH -D ./
#SBATCH -o ./results/perf/%j.out
#SBATCH -e ./results/perf/%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
 
module load slurm_setup

perf record --call-graph=dwarf -F 997 -o perf.data ./MolSim-WS25-GroupD/build/project/MolSim ./input/benchmark_10k.yaml
