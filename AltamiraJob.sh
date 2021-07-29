#!/bin/bash
#
#SBATCH --job-name=eclipse
#SBATCH --output=eclipse.out
#SBATCH --error=eclipse.err
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=3500
#SBATCH --time=00:05:00
# From here the job starts
mpirun ./ECLIPSE_TEB 2 Example_Folder/Example.ini
