suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(genomeIntervals)) ## Used for readGff3
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(argparser))

#### Read command-line arguments ####
argp <- arg_parser(paste(
    "Use a reference genomic-features file (.gff)",
    "to find genes in a region of interest"
))
argp <- add_argument(argp, "--samplegroup",
                     help = "Group name of related samples"
)
argp <- add_argument(argp, "--refpath",
                     help = "file path of reference"
)
argp <- add_argument(argp, "--refstrain",
                     help = "strain name for reference"
)
argp <- add_argument(argp, "--region",
                     type = "character",
                     help = "Regions are specified as: seqname[:STARTPOS[-ENDPOS]] e.g. 08:409282-429234"
)

argv <- parse_args(argp)

refDir <- argv$refpath
pfCurrRefs <- data.frame(
    strain = c("3D7", "Dd2", "Supp"),
    version = c("52", "57", "52")
)
ref <- filter(pfCurrRefs, strain == argv$refstrain)
## Not looking for regions of interest in supplementary sequences yet

# ---------- Read genomic features --------------------------------------------
file.gff <- file.path(
    refDir,
    paste0(
        "PlasmoDB-", ref$version, "_Pfalciparum", ref$strain,
        ".gff"
    )
)

pf_features <- tryCatch(
    readGff3(file.gff, quiet = TRUE),
    error = function(e) {
        stop(paste(e, "Filepath:", file.gff))
    }
)
# ---------- Find features in region of interest ------------------------------
chrname <- str_remove(argv$region, ":.+")
gstart <- str_remove(argv$region, ".+:") |>
    str_remove("-.+")
gend <- str_remove(argv$region, ".+:") |>
    str_remove(".+-")
roigr <- GRanges(
    seqnames = levels(seqnames(pf_features)) |>
        str_subset(paste0("Pf", ref$strain, "_", chrname)),
    ranges = IRanges(start = gstart, end = gend)
)
featROIidx <- findOverlaps(GRanges(pf_features), roigr,
                           type = "within", select = "last"
)
featROI <- pf_features[!is.na(featROIidx)]
# ---------- Filter to genes and write to csv ------------------------------
geneROI <- featROI[
    annotation(featROI)$type == "protein_coding_gene",
]

data.frame(ID = getGffAttribute(geneROI, "ID"),
           GeneName = getGffAttribute(geneROI, "Name"),
           Description = getGffAttribute(geneROI, "description") |>
               str_replace_all("\\+", " ")
) |>
    write_csv(paste0(argv$samplegroup, "genes_chr", chrname, "_ROI.csv"))
