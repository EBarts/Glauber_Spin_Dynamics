#!/bin/bash
#SBATCH --time=28:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=FSweep20
#SBATCH --output=CPPtestPar.out
#SBATCH --mem=64000

module load foss/2020a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./FsweepParOMP