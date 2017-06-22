#!/bin/bash
#
#SBATCH --job-name=ASFracture
#SBATCH --output=ASFracture_output.txt
#SBATCH --error=ASFracture_error.txt
#SBATCH --nodelist=aquinonode1  --ntasks=64
#SBATCH --time=4800:00:00





mpiexec ../../ASFracture-opt -i asfractureMaster.i
