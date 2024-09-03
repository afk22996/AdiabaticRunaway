#!/bin/bash
#SBATCH --job-name penguin
#SBATCH --nodes 1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 128gb
#SBATCH --time 72:00:00
#SBATCH --gpus-per-node a100:2

module load cuda
make clean
make
./penguin