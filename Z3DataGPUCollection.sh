#!/bin/bash

#SBATCH -p gpu --gres=gpu:4

#SBATCH -n 2
#SBATCH -t 30:00:00

module load cuda

./z3main
