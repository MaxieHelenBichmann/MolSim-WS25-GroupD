#!/bin/bash
#SBATCH -J intel-omp-profile
#SBATCH -D ./
#SBATCH -o ./results/intel/%j.out
#SBATCH -e ./results/intel/%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00

# Load necessary modules for profiling
module load slurm_setup
module load intel-oneapi-vtune/2024.2.1 
module load intel-oneapi-advisor/2024.2.1
module load gcc/15.2.0

# Set OpenMP threads to match the allocated CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"

mkdir -p results/intel
make -C MolSim-WS25-GroupD/build vtune-hotspots
tar -czf MolSim-WS25-GroupD/build/profiling-results.tar.gz -C MolSim-WS25-GroupD/build vtune/hotspots
