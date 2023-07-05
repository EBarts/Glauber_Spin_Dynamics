#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=FSweepSpinSave
#SBATCH --output=CPPtestPar.out
#SBATCH --mem=64000

module load foss/2020a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./FsweepParOMP