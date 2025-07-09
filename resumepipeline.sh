#!/bin/bash
echo Additional options can be added to nextflow by enclosing in quotes.
echo E.g. ./resumepipeline \'sad_nobel -c path-to-task-config\'

module purge
module load nextflow/25.04.2
nextflow run main.nf -profile apptainer -resume $1
