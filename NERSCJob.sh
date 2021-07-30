#!/bin/bash
#SBATCH --job-name=eclipse
#SBATCH --output=eclipse.out
#SBATCH --error=eclipse.err
#SBATCH --qos=regular
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --constraint=haswell
#SBATCH --mail-user=<user>@<domain>.com
#SBATCH --mail-type=ALL

srun ./ECLIPSE_TEB 2 Example/Example.ini
