#!/bin/bash
#SBATCH --job-name closestFlows
#SBATCH --nodes 1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 80gb
#SBATCH --time 72:00:00

module load anaconda3
python closestFlow.py