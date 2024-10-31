## Read gridss vcf that has been filtered with gridss_somatic_filter.R
## and filter by requiring that a majority of resistant samples pass quality
## criteria.
## Alternative filtering requires a majority of resistant samples to have AF
## greater than some critical value, and the maximum AF among resistant samples
## to be greater than the maximum among parent samples

#### required libraries ####
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library(VariantAnnotation))
library(argparser)

#### Read command-line arguments ####
argp <- arg_parser(paste(
  "Read an SV vcf filtered to 'somatic' events, and filter again."
))
argp <- add_argument(argp, "--samplegroup",
  help = "Group name of related samples. Required"
)
argp <- add_argument(argp, "--scriptdir",
  help = "Path to libgridss.R script"
)

argp <- add_argument(argp, "--parentlist",
  help = "parent filenames in a space separated string"
)
argp <- add_argument(argp, "--critsamplecount",
  type = "double",
  help = paste(
    "Minimum number of (non-parent) samples",
    "having an SV. Default is (#samples + 1)/2"
  )
)
argv <- parse_args(argp)


## GRIDSS function readVcf copied from "libgridss.R":
readVcf <- function(file, ...) {
  raw_vcf <- VariantAnnotation::readVcf(file = file, ...)
  if (!all(unlist(alt(raw_vcf)) != "")) {
    write(
      "Performing work-around for https://github.com/Bioconductor/VariantAnnotation/issues/8",
      stderr()
    )
    alt <- read_tsv(
      file,
      comment = "#",
      col_names = c(
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "INFO", "FORMAT", seq_len(ncol(geno(raw_vcf)[[1]]))
      ),
      cols_only(ALT = col_character())
    )$ALT
    VariantAnnotation::fixed(raw_vcf)$ALT <- CharacterList(lapply(as.character(alt), function(x) x))
  }
  raw_vcf <- fix_parid(raw_vcf)
  if ("MATEID" %in% names(info(raw_vcf)) & !is.null(info(raw_vcf)$MATEID)) {
    raw_vcf <- flatten_mateid(raw_vcf)
  }
  return(raw_vcf)
}
flatten_mateid <- function(vcf) {
  if (!("MATEID" %in% names(info(vcf)))) {
    stop("Missing MATEID")
  }
  mateid <- info(vcf)$MATEID
  if (any(elementNROWS(mateid)) > 1) {
    stop("Multiple MATEID for a single record not supported.")
  }
  info(vcf)$MATEID <- as.character(mateid)
  return(vcf)
}
fix_parid <- function(vcf) {
  if (
    ("PARID" %in% row.names(info(header(vcf))) & !("MATEID" %in% row.names(info(header(vcf))))) |
      (!is.null(info(vcf)$PARID) & any(elementNROWS(info(vcf)$PARID) > 0) &
        (is.null(info(vcf)$MATEID) | all(elementNROWS(info(vcf)$MATEID) == 0)))) {
    parid <- info(vcf)$PARID
    info(vcf)$PARID <- NULL
    parid_header_ordinal <- which(row.names(info(header(vcf))) == "PARID")
    row.names(info(header(vcf)))[parid_header_ordinal] <- "MATEID"
    info(header(vcf))$Description[parid_header_ordinal] <- "ID of mate breakends"
    info(vcf)$MATEID <- parid
    write("WARNING: MATEID header not found. Assuming VCF was generated prior to GRIDSS 2.8.0 and rewriting PARID as MATEID.", stderr())
  }
  return(vcf)
}

#### Read input vcf ####
somvcfname <- paste0(
  argv$samplegroup,
  "_high_and_low_confidence_somatic.vcf.bgz"
)
somvcf <- tryCatch(
  readVcf(somvcfname),
  error = function(e) {
    stop(paste(
      e, "when reading GRIDSS somatic vcf", somvcfname
    ))
  }
)
#----------- Filtering parameters -----------------
parentlist <- str_split(argv$parentlist, " ") |>
  unlist() |>
  str_remove("_nodup.bam$")
samplesOI <- setdiff(rownames(colData(somvcf)), parentlist)
## For replicate samples, a filtering criterion is that variants are
## present in critSamplesSom or more of the samples.
## The default value is half the number of (non-parent) samples.
if (is.na(argv$critsamplecount)) {
  critSamplesSom <- (1 + length(samplesOI)) / 2
} else {
  critSamplesSom <- argv$critsamplecount
}

BQcrit <- 200
Qcrit <- 200
acceptableFilt <- c("PASS", "imprecise")
firstFiltname <- paste0(argv$samplegroup, "_high_and_imprecise_somatic.vcf")

AFcrit <- 0.15

#### Filter by event has FILTER values in acceptable set ####
passvcf <- VariantAnnotation::subset(somvcf, filt(somvcf) %in% acceptableFilt)
## include mates
passvcf <- VariantAnnotation::subset(
  somvcf, names(somvcf) %in% c(names(passvcf), na.omit(info(passvcf)$MATEID))
)
rm(somvcf)
writeVcf(passvcf, firstFiltname)

#### Filter by event is in 'majority' of samples and not parent ####
majsom <- passvcf[
  apply(
    geno(passvcf)$BQ, 1,
    function(Q) {
      sum(Q[samplesOI] > BQcrit) >= critSamplesSom
    }
  ) &
    apply(
      geno(passvcf)$QUAL, 1,
      function(Q) {
        sum(Q[samplesOI] > Qcrit) >= critSamplesSom
      }
    )
]

somMinAF <- passvcf[
  apply(
    geno(passvcf)$AF, 1,
    function(af) {
      sum(af[samplesOI] |> unlist() > AFcrit) >= critSamplesSom &
        max(af[samplesOI] |> unlist()) > 2 * max(unlist(af[parentlist]))
    }
  )
]

#### Save the union of the 2 filters as vcf and summary table ####
somEitherFilt <- rbind(
    majsom, somMinAF
) |> unique() |> sort()
writeVcf(somEitherFilt, 
         paste0(argv$samplegroup, "_somatic_by_QUALorAF.vcf"))

filtdf <- data.frame(
  gridssID = rownames(somEitherFilt),
  seqname = seqnames(somEitherFilt),
  pos = start(somEitherFilt),
  REF = ref(somEitherFilt) |> unlist(),
  ALT = alt(somEitherFilt) |> unlist() |> str_trunc(width = 24, side = "right"),
  as.data.frame(geno(somEitherFilt)$AF) |> unnest(cols = everything()) |>
    rename_with(~ paste0("AF_", .x)),
  as.data.frame(geno(somEitherFilt)$QUAL) |> rename_with(~ paste0("QUAL_", .x))
) |>
    arrange(gridssID)
write_csv(
  filtdf,
  file.path(paste0(argv$samplegroup, "_somatic_by_QUALorAF.csv"))
)
