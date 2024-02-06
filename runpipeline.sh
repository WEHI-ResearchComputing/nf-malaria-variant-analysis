#!/bin/bash
module purge
module load miniconda3/latest nextflow/23.10.0  apptainer/1.0.0 
nextflow run main.nf -profile milton  -config config/nextflow_sh.config 
