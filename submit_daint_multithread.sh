#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --account=pr133
#SBATCH --time=00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=multithread
#SBATCH --output=test_fast.o
#SBATCH --error=test_fast.e

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -u python precipitation_indices.py
