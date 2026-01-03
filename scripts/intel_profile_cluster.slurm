#!/bin/bash
#SBATCH -J=intel-profile
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


module load intel-oneapi-vtune/2024.2.1 
module load intel-oneapi-advisor/2024.2.1

make intel-profile-all

# Get final archive path
ARCHIVE_PATH="$(pwd)/profiling-results.tar.gz"
ARCHIVE_SIZE=$(du -h $ARCHIVE_PATH | cut -f1)

echo ""
echo "=========================================="
echo "Profiling Complete!"
echo "=========================================="
echo "Results archive: $ARCHIVE_PATH"
echo "Archive size: $ARCHIVE_SIZE"
echo "Job ID: $SLURM_JOB_ID"
echo "End time: $(date)"
echo ""
echo "To download to your local machine:"
echo "  scp $(hostname):$ARCHIVE_PATH ."
echo ""
echo "To view results locally:"
echo "  tar xzf profiling-results.tar.gz"
echo "  vtune-gui profiling-results/vtune/hotspots"
echo "  advisor-gui profiling-results/advisor/survey"
echo "=========================================="
