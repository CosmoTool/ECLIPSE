#!/bin/bash
#SBATCH --job-name=ECLIPSE
#SBATCH --output=ECLIPSE64.out
#SBATCH --error=ECLIPSE64.err
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=3
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH --mail-user=bilbao@ifca.unican.es
#SBATCH --mail-type=ALL

srun ./ECLIPSE_TEB 12 NSide64/TEB.ini
