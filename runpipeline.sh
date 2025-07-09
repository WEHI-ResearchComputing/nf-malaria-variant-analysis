#!/bin/bash
module purge
echo Additional options can be added to nextflow by enclosing in quotes.
echo E.g. ./runpipeline \'-c path-to-task-config\'
module load nextflow/25.04.2
nextflow run main.nf -profile apptainer $1
