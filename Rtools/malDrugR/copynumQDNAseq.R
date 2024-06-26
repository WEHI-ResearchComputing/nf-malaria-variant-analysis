## This reads the reference genome of P. falciparum in BioStrings format,
## reads bam files (first time) or binned read counts if available,
## then uses functions from QDNAseq to calculate read depths and copynumbers.
## A 'mappability file' is used for correcting the counts
##
#Rscript ${projectDir}/bin/malDrugR/copynumQDNAseq.R \
    #    --samplegroup ${groupId} \
    #    --samplekeyfile ${samplekeyfile} \
    #    --groupkeyfile ${groupkeyfile} \
    #    --bin_in_kbases ${params.bin_in_kbases}
#### Read command-line arguments ####
library(argparser)
argp <- arg_parser(paste(
  "Read a set of bam files, ", "count reads in bins defined on a reference genome",
  ", and calculate copy number scores. Intermediate files are read if available.",
  "\nData frames of normalised CN and CN of samples scaled relative to a 'parent'",
  " sample are saved as RDS."
))
argp <- add_argument(argp, "--samplegroup",
                     help = "Group name of related samples. Required"
)
argp <- add_argument(argp, "--parentId",
                     help = "Parent ID of the group."
)
argp <- add_argument(argp, "--bams",
                     help = "file containing group info including short ID for parent."
)
argp <- add_argument(argp, "--bin_in_kbases",
                     default = "1",
                     help = "bin size in kbp"
)
argp <- add_argument(argp, "--refDir",
                     help = "reference directory "
)
argp <- add_argument(argp, "--bsref",
                     help = "BS reference name to import library"
)
argv <- parse_args(argp)


#### required libraries ####
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("QDNAseq"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library(devtools))

groupId <- argv$samplegroup
parentID <- argv$parentId
refDir <- argv$refDir
groupbamL <- strsplit(argv$bams, " ")[[1]]
groupbamL <- groupbamL[!grepl(parentID, groupbamL)]
sampleL <- sub("_nodup\\.bam$", "", groupbamL)

library(argv$bsref,
  character.only = TRUE
)
pfg <- get(argv$bsref)


#### Function to create an annotated dataframe of bins, ####
##   listing their start, end, GC percentage, etc.
##   Exclude non-nuclear
makePfBins <- function(bin_in_kbases) {
  binRfile <- file.path(refDir, paste0("PfBins", bin_in_kbases, "kb.rds"))
  if (file.exists(binRfile)) {
    pfBins <- readRDS(binRfile)
  } else {
    pfBins <- createBins(pfg, as.numeric(bin_in_kbases)
    )
    mapityfile <- file.path(
      refDir, "mappability", "kmer30_err2.bw"
    )
    if( !file.exists(mapityfile) ){
      stop(paste("Mappability file", mapityfile, "not found"))}
    ### Untested! bigWig executables require module ucsc-tools
    remotes::install_github("WEHI-ResearchComputing/EnvironmentModules")
    library(EnvironmentModules)
    module_load("ucsc-tools/331")
    bigWigpath <- file.path(
      "/stornext/System/data/apps/ucsc-tools/ucsc-tools-331/bin",
      "bigWigAverageOverBed"
    )
    if( !file.exists(bigWigpath) ){
      stop(paste(bigWigpath, "not found."))}
    pfBins$mappability <- calculateMappability(
      pfBins,
      bigWigFile = mapityfile,
      bigWigAverageOverBed = bigWigpath,
      chrPrefix = ""
    )
    pfBins <- AnnotatedDataFrame(
      pfBins,
      varMetadata = data.frame(labelDescription = c(
        "Chromosome name",
        "Base pair start position",
        "Base pair end position",
        "Percentage of non-N nucleotides (of full bin size)",
        "Percentage of C and G nucleotides (of non-N nucleotides)",
        "Average mappability of 30mers"
      ), row.names = colnames(pfBins))
    )
    saveRDS(pfBins, binRfile)
  }
}
###### Binned counts for samples ######

## Parent sample handled separately because many groups use the same parent
parentcountsn <- file.path(
  paste0("Counts", argv$bin_in_kbases, "k_", parentID, ".rds")
)

if (file.exists(parentcountsn)) {
  counts_parents <- readRDS(parentcountsn)
} else {
  if (!exists(deparse(substitute(pfBins)))) {
    pfBins <- makePfBins(argv$bin_in_kbases)
  }
  counts_parents <- binReadCounts(
    bins = pfBins,
    bamfiles = c(file.path(paste0(parentID, ".bam"))),
    bamnames = parentID
  )
  saveRDS(counts_parents, parentcountsn)
}


## Resistant samples
countfilen <- paste0(
  "Counts", argv$bin_in_kbases, "k_", groupId,
  ".rds"
)
if (file.exists(file.path(countfilen))) {
  countsinbins <- readRDS(file.path(countfilen))
} else {
  if (!exists(deparse(substitute(pfBins)))) {
    pfBins <- makePfBins(argv$bin_in_kbases)
  }
  countsinbins <- binReadCounts(
    pfBins,
    bamfiles = file.path(groupbamL),
    bamnames = sampleL
  )

  saveRDS(
    countsinbins,
    file.path(countfilen)
  )
}

countsinbins <- Biobase::combine(counts_parents, countsinbins)

## Filtering with min-mappability=50% and NO blacklist
## No residuals calculated
countsFiltered <- applyFilters(
  countsinbins,
  residual = FALSE, mappability = 50, blacklist = FALSE
)

## Estimate the correction for GC content and mappability.
## This calculates a loess fit to the counts
## (as a function of GC and mappability) and adds it to the data structure
countsCorrected <- estimateCorrection(countsFiltered)

## Use the values added by estimateCorrection to correct the data:
## Copy numbers are counts divided by 'fit'
copyNums <- correctBins(countsCorrected)
## Scale copy numbers by median (default)
copyNumsNormed <- normalizeBins(copyNums)
## Scale strain copy numbers by dividing by parent copy numbers
StrainScaledByParents <- compareToReference(
  copyNumsNormed,
  c(FALSE, rep(1, times = length(groupbamL)))
)
## Function to convert QDNASeq object to data frame
convertQDNAtoDF <- function(qobject) {
  if ("counts" %in% assayDataElementNames(qobject)) {
    scoreDF <- as.data.frame(assayDataElement(qobject, "counts"))
  } else if ("copynumber" %in% assayDataElementNames(qobject)) {
    scoreDF <- as.data.frame(assayDataElement(qobject, "copynumber"))
  }
  scoreDF$chrom <- sapply(strsplit(rownames(scoreDF), ":"),
                          function(x) (unlist(x)[1]))
  scoreDF$range <- sapply(strsplit(rownames(scoreDF), ":"),
                          function(x) (unlist(x)[2]))
  scoreDF$start <- sapply(strsplit(scoreDF$range, "-"),
                          function(x) (as.numeric(unlist(x)[1])))
  scoreDF$end <- sapply(strsplit(scoreDF$range, "-"),
                        function(x) (as.numeric(unlist(x)[2])))

  return(scoreDF)
}
CN_df <- convertQDNAtoDF(copyNumsNormed)
saveRDS(
  CN_df,
  file.path(paste0(groupId, ".CN_df_",
      argv$bin_in_kbases, "k", ".rds"
    )
  )
)
scaled_df <- convertQDNAtoDF(StrainScaledByParents)
saveRDS(
  scaled_df,
  file.path(paste0(groupId, ".CN_compare_df_",
      argv$bin_in_kbases, "k", ".rds"
    )
  )
)
