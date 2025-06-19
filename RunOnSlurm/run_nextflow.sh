#!/bin/bash
#SBATCH --job-name=nextflow_master
#SBATCH --output=nextflow_master_%j.log
#SBATCH --error=nextflow_master_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short

# Load necessary modules
module load nextflow
module load singularity
module load java
module load wget

# Set Java home dynamically using which java
export JAVA_HOME=$(dirname $(dirname $(which java)))
export PATH=$JAVA_HOME/bin:$PATH

# Set Nextflow Java home
export NXF_JAVA_HOME=$JAVA_HOME

# Run Nextflow
chmod +x bin/*
nextflow run main.nf -profile slurm  --input samplesheet.tsv --outdir ./outputs/