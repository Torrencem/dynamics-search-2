#!/bin/bash

#SBATCH -p gpu --gres=gpu:8

#SBATCH -n 2
#SBATCH -t 30:00:00

#SBATCH -e e-z4.out

module load cuda

./build.sh

./main
