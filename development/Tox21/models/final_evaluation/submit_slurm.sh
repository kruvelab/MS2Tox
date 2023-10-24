#!/bin/bash

#SBATCH --time=120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH --job-name="Tox21_models"
#SBATCH --partition=amd

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.5.3
module load squashfs/4.4

nextflow run test_models.nf