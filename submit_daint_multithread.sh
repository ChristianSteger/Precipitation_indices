#!/bin/bash -l
#SBATCH --job-name=precipitation_indices
#SBATCH --account=pr133
#SBATCH --time=00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=multithread
#SBATCH --output=precipitation_indices.o
#SBATCH --error=precipitation_indices.e

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -u python precipitation_indices.py
