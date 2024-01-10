# nf-malaria-variant-analysis

## To run on Milton
```
module load anaconda3 nextflow/23.10.0  apptainer/1.0.0 
nextflow run main.nf -profile milton  -config config/nextflow_sh.config
```

This will use data files from `/vast/projects/malaria/data_sample/sh_16000` which are samples from the original data files.

