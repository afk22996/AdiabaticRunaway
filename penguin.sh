#/bin/bash
#PBS -N penguin
#PBS -l select=1:ncpus=1:mem=20gb:gpu_model=v100:ngpus=2
#PBS -l walltime=72:00:00
#PBS -o output.log
#PBS -j oe

cd $PBS_O_WORKDIR
module load cuda/12.0.1-gcc
./penguin