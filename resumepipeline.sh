#!/bin/bash
module purge
module load nextflow/23.10.0
nextflow run main.nf -profile apptainer -resume
