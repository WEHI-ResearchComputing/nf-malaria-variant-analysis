# nf-malaria-variant-analysis

## To run on Milton
```
./runpipeline.sh
```
## To resume on Milton
```
./resumepipeline.sh
```

This will use all config parameters found in `nextflow.config`.

## Format of samplegroup.txt

It is a tab-delimited text file, with 5 columns. The columns are (`groupId`, `sampleId`, `fastqbase`, `ref`, `parentId`)

### **Notes**
* `sampleid` must contain `groupid`
* No need for parentbamlist
* If parent is not a sample add a line with no parentid

