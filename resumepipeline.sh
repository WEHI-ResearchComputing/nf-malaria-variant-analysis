#!/bin/bash
module purge
module load anaconda3 nextflow/23.10.0  apptainer/1.0.0 
nextflow run main.nf -profile milton  -config config/nextflow_sh.config -resume