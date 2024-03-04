## Reading vcf files output by bcftools mpileup+call variant caller
## Filter to new-in-resistant-strain, and in-gene and reliable

# ---------- Read arguments ---------------------------------------------------
#
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(argparser))
#### Load remaining libraries ####
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(genomeIntervals)) ## Used for readGff3
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(AnnotationDbi))

#### Read command-line arguments ####
argp <- arg_parser(paste(
  "Read a small-variant vcf file written by bcftools",
  "and filter to good-quality variants that are not in reference samples"
))
argp <- add_argument(argp, "--samplegroup",
  help = "Group name of related samples. Required"
)
argp <- add_argument(argp, "--refpath",
  help = "file path of reference"
)
argp <- add_argument(argp, "--refstrain",
  help = "strain name for reference"
)
argp <- add_argument(argp, "--parentlist",
  help = "parent filenames in a space separated string"
)

argp <- add_argument(argp, "--QUALcrit",
  default = 50,
  help = paste(
    "Critical value for filtering variants by",
    "QUAL of variant call"
  )
)
argp <- add_argument(argp, "--critsamplecount",
  type = "double",
  help = paste(
    "Minimum number of non-parent samples",
    "having a variant. Default is (#samples + 1)/2"
  )
)
argv <- parse_args(argp)



#### Read group details and set file paths ####

parentlist <- sub("_nodup\\.bam$", "", strsplit(argv$parentlist, " ")[[1]])
refDir <- argv$refpath

varDir <- "./"

pfCurrRefs <- data.frame(
  strain = c("3D7", "Dd2"),
  version = c("52", "57")
)
ref <- filter(pfCurrRefs, strain == argv$refstrain)

# ---------- Read genomic features --------------------------------------------
file.gff <- file.path(
  refDir,
  paste0("PlasmoDB-", ref$version, "_Pfalciparum", ref$strain, ".gff")
)
pf_features <- tryCatch(
  readGff3(file.gff, quiet = TRUE),
  error = function(e) {
    stop(paste(e, "Filepath:", file.gff))
  }
)
## Filter to remove PfEMP1 var, rifin, and STEVOR genes, and pseudogenes
## - genes known to be highly variable and not relevant
idxGene2remove <- grepl(
  "PfEMP1|rifin|stevor|pseudogene",
  annotation(pf_features)$gffAttributes
)
nameGene2remove <- pf_features[idxGene2remove]$gffAttributes %>%
  str_extract(., paste0("P[Ff]", ref$strain, "_[0-9]{7}")) %>%
  unique()
## None of these genes are on mitochondria (PfDd2_MT.* or PF3D7_MIT.*)
idx2remove <- str_detect(
  annotation(pf_features)$gffAttributes,
  paste(nameGene2remove, collapse = "|")
)
pf_featuresNovar <- pf_features[!idx2remove]
CDSnoVar <- pf_featuresNovar[annotation(pf_featuresNovar)$type == "CDS", ]
file.chrinfo <- file.path(
  refDir, paste0(
    "PlasmoDB-", ref$version, "_Pfalciparum", ref$strain,
    "_Genome.fasta.fai"
  )
)
chrinfo <-
  tryCatch(
    read.table(
      file = file.chrinfo,
      sep = "\t",
      col.names = c(
        "seqnames", "seqlengths",
        "OFFSET", "LINEBASES", "LINEWIDTH"
      )
    ),
    error = function(e) {
      stop(paste(e, "Filepath:", file.chrinfo))
    }
  )
grPfalc <- GRanges(
  seqnames = chrinfo$seqnames,
  ranges = IRanges(start = 1, width = chrinfo$seqlengths)
)

if (ref$strain == "3D7") { # haven't yet got this data for Dd2
  corepath <- file.path(refDir, "GenomeRegionsMilesEtAL_plus_10centromere.txt")
  getRegions <- function(filep) {
    coreGenome <- read_delim(
      filep,
      delim = " ",
      col_types = "ciici"
    ) %>% filter(Type == "Core")
    grCore <- GRanges(
      seqnames = coreGenome$Chromosome,
      ranges = IRanges(
        start = coreGenome$Start,
        end = coreGenome$Stop,
        width = coreGenome$Size
      ) # specify all 3 as sanity check against bed-format confusion
    )
    ## Extend 'core genome' to include mitochondrial
    return(c(grCore, grPfalc[which(seqnames(grPfalc) == "Pf3D7_MIT_v3")]))
  }
  grCoreMT <- tryCatch(
    getRegions(corepath),
    error = function(err) {
      warning(paste(
        err, corepath,
        "\nRevert to full genome instead of core regions."
      ))
      return(grPfalc)
    }
  )
}

# ---------- Define filter functions -------------------------------------------
#
filt_vcf <- function(vcf, QUALcrit,
                     parents = "parent") {
  ## Discriminate between parent samples and samples of interest
  samplesOI <- setdiff(rownames(colData(vcf)), parents)
  ## Filter by: overall QUAL is above threshold
  vcf_qualfilt <- vcf[
    rowRanges(vcf)$QUAL > QUALcrit
  ]
  ## Filter by: not an InDel of AT run or repetitive region.
  vcf_ATfilt <- vcf_qualfilt[
    !info(vcf_qualfilt)$INDEL |
      sapply(
        names(rowRanges(vcf_qualfilt)),
        function(n) {
          vindel <-
            str_remove(
              n, paste0("P[fF]", ref$strain, ".*[0-9]+_")
            ) %>%
            str_split_fixed(., "/", 2)
          changeset <- union(
            # remove 1st 'anchor' nucleotide
            vindel[1] %>%
              str_sub(2, -1) %>%
              str_extract_all(., "") %>%
              unlist(),
            vindel[2] %>% str_sub(2, -1) %>%
              str_extract_all(., "") %>% unlist()
          )
          allAT <- !"G" %in% changeset & !"C" %in% changeset
          # indel is contained in ref, or vice versa
          containsAlt <- str_detect(vindel[1], vindel[2]) |
            str_detect(vindel[2], vindel[1])
          containsRun <- any(str_detect(
            vindel, "AAAAAA|TTTTTT|ATATAT"
          ))
          ## Combined filter rule: keep only if short or complex
          max(nchar(vindel)) < 5 |
            (!allAT & !containsRun & !containsAlt
            )
        }
      )
  ]
  ## Filter by: there are genotype calls in the samples of interest
  ## that are not in the set of parent genotypes plus reference (GT 0).
  vcf_somatic <- vcf_ATfilt[
    apply(geno(vcf_ATfilt)$GT, 1, function(GT) {
      setdiff(GT[samplesOI], c("0", ".", GT[parents])) %>%
        length(.) > 0
    })
  ]

  ## Filter by: event is in core genome or mitochondria
  if (ref$strain == "3D7") { # haven't yet got this data for Dd2
    vcf_posfilt <- vcf_somatic[
      overlapsAny(vcf_somatic, grCoreMT)
    ]
    return(vcf_posfilt)
  } else {
    return(vcf_somatic)
  }
}

## Function to get list of new 'somatic' Alt IDs
somatic <- function(gt, pgt) {
  gtdf <- as.data.frame(gt)
  filter(gtdf, !gt %in% c("0", ".", pgt))
}

# ---------- Read vcf files, and filter --------------------------------
#
samplevcf <- tryCatch(
  readVcf(paste0(argv$samplegroup, ".vcf")),
  error = function(e) {
    stop(paste(e, "Filepath:", paste0(argv$samplegroup, ".vcf")))
  }
)
## Filter
eventFilt <- filt_vcf(
  vcf = samplevcf, parents = parentlist,
  QUALcrit = argv$QUALcrit
)

## Filter again to majority of variants in samples are absent in parents
samplesOI <- setdiff(rownames(colData(samplevcf)), parentlist) |>
  sort()
## For replicate samples, a filtering criterion is that variants are
## present in critSamplesSom or more of the samples.
## The default value is half the number of (non-parent) samples.
if (is.na(argv$critsamplecount) || length(argv$critsamplecount) == 0) {
  critSamplesSom <- (1 + length(samplesOI)) / 2
} else {
  critSamplesSom <- argv$critsamplecount
}

majsom <- eventFilt[
  apply(
    geno(eventFilt)$GT, 1,
    function(GT) {
      nrow(somatic(GT[samplesOI], GT[parentlist])) >= critSamplesSom
    }
  )
]
writeVcf(
  majsom, file.path(
    varDir,
    paste0(argv$samplegroup, critSamplesSom, "plus", ".Qcrit", argv$QUALcrit, ".vcf")
  )
)

#### Events that are in gene features in GFF ####
## indels in introns and UTRs are potentially of interest.
## Single nucleotide events are only interesting in CDS

#### Non-synonymous SNPs in coding regions  ####
eventCDS <- findOverlaps(rowRanges(majsom), GRanges(CDSnoVar))
snpCDS <- subset(majsom[queryHits(eventCDS)], !INDEL)

#### Read bam files to get allele counts of snpCDS positions ###
countalleles <- function(bamfile, vcf) {
  bamparams <- ScanBamParam(
    which = rowRanges(vcf),
    what = c("qname", "pos", "seq", "qual")
  )
  Reads <- readGAlignments(
    file = bamfile,
    param = bamparams,
    use.names = TRUE
  )
  readCounts <- map(names(vcf), function(Name) {
    vcfgr <- rowRanges(vcf)[Name]
    width1stAlt <- unlist(vcfgr$ALT)[1] |> width()
    eventSeq <- Reads[
      findOverlaps(vcfgr, Reads, type = "within") |>
        subjectHits()
    ]
    PosInSeq <- mapToAlignments(vcfgr, eventSeq)
    refwidthDNA <- 
      subseq(mcols(eventSeq)$seq,
        start = start(PosInSeq),
        width = width(vcfgr$REF)
    )
    altwidthDNA <- tryCatch(
      subseq(mcols(eventSeq)$seq,
        start = start(PosInSeq),
        width = width1stAlt
      ),
      error = function(e) {
        print(paste(vcfgr, "ALT does not fit within reads. Skipping"))
      }
    )
    countresult <- data.frame(
      SampleName = character(), variant = character(),
      TotCount = integer(), RefCount = integer(),
      AltCount = integer()
    )
    if (class(altwidthDNA) == "DNAStringSet") {
      countresult <- data.frame(
        SampleName = basename(bamfile) |>
          str_remove("_s.bam") |> str_remove("_nodup.bam"),
        variant = Name,
        TotCount = length(eventSeq),
        RefCount = sum(refwidthDNA == vcfgr$REF),
        AltCount = sum(altwidthDNA == unlist(vcfgr$ALT)[1])
      )
    }
    return(countresult)
  }) |> list_rbind()
}

SNPalleleCounts <- map(
  c(
    paste0(samplesOI, "_nodup.bam"),
    paste0(argv$parentlist)
  ),
  countalleles,
  vcf = snpCDS
) |>
  list_rbind()


SNPalleleCounts <- arrange(SNPalleleCounts, variant)

#### Annotate SNPs ####
#### Load genome and transcript db
transcriptdb <- file.path(
  refDir,
  paste0("PlasmoDB-", ref$version, "_Pfalciparum", ref$strain, "_txdb.sql")
)
txdb <- tryCatch(
  txdb <- loadDb(transcriptdb),
  error = function(e) {
    stop(paste(e, "Filepath:", transcriptdb))
  }
)
library(paste0("BSgenome.Pfalciparum", ref$strain, ".PlasmoDB.", ref$version),
  character.only = TRUE
)
pfg <- get(paste0("BSgenome.Pfalciparum", ref$strain, ".PlasmoDB.", ref$version))

#### AA prediction for SNPs in CDS
AApred <- predictCoding(query = snpCDS, subject = txdb, seqSource = pfg)

#### SNPs in same codon
mcols(AApred)$warning <- NA
if (paste(AApred$CDSID, AApred$PROTEINLOC) |> unlist() |> n_distinct() <
  length(AApred)
) {
  multihit <- which(mcols(AApred) |>
    as.data.frame() |>
    add_count(CDSID, PROTEINLOC) |>
    dplyr::select(n) > 1)
  mcols(AApred[multihit])$warning <- "multihit on codon"
}

#### Remove synonymous ####
AApred <- AApred[which(mcols(AApred)$CONSEQUENCE != "synonymous" |
  !is.na(mcols(AApred)$warning))]
nonsynSNP <-
  cbind(
    as.data.frame(AApred) |> dplyr::select(-ALT),
    data.frame(
      SNP = names(AApred),
      ALT = as.character(unlist(AApred$ALT)),
      AAchanges = with(
        AApred,
        paste0(REFAA, PROTEINLOC, VARAA)
      )
    )
  ) |> dplyr::select(SNP, seqnames,
    pos = start, strand, REF, ALT, QUAL,
    AAchanges, REFCODON, VARCODON, warning
  )

if (all(is.na(nonsynSNP$warning))) {
  nonsynSNP <- dplyr::select(nonsynSNP, !c(REFCODON, VARCODON, warning))
}
SNPdetails <- left_join(
  nonsynSNP,
  SNPalleleCounts |> mutate(
    AltFrac = paste0(AltCount, "/", TotCount),
    TotCount = NULL, RefCount = NULL, AltCount = NULL
  ) |>
    pivot_wider(
      names_from = c(SampleName),
      values_from = c(AltFrac),
      names_prefix = "AF_"
    ),
  by = join_by("SNP" == "variant")
) |>
  dplyr::select(!SNP)

#### Get gene details for indels and non-synonymous SNPs ####
eventGene <- findOverlaps(rowRanges(majsom), GRanges(pf_featuresNovar),
  select = "last"
)
indelGene <- subset(majsom[!is.na(eventGene)], INDEL)
if (nrow(indelGene) > 0) {
  indelGeneAttr <- pf_featuresNovar[
    findOverlaps(rowRanges(indelGene), GRanges(pf_featuresNovar),
      select = "last"
    )
  ]
  #### Get gene details for indels ####
  indels.Feat.df <- data.frame(
    rowRanges(indelGene)[, c("REF", "ALT", "QUAL")],
    ID = getGffAttribute(indelGeneAttr, "ID") |> unlist()
  )

  baseIDEvents <- str_replace(
    indels.Feat.df$ID, "^.*P", "P"
  ) |>
    str_remove("[\\.\\-].*$")
  geneDetail <- map(baseIDEvents, function(ID) {
    genegff <- pf_featuresNovar[getGffAttribute(pf_featuresNovar, "ID") == ID]
    data.frame(
      GeneName = getGffAttribute(genegff, "Name"),
      Description = getGffAttribute(genegff, "description")
    )
  }) |>
    list_rbind()
  indels.Feat.df <- cbind(indels.Feat.df, geneDetail) |>
    dplyr::select(-width, -strand) |>
    mutate(
      Gene = case_when(
        is.na(GeneName) ~ Description |>
          str_replace_all(pattern = "\\+", replacement = " "),
        TRUE ~ GeneName
      ),
      ALT = CharacterList(ALT) |> unstrsplit(sep = ","),
      GeneName = NULL, Description = NULL
    )
} else {
  indels.Feat.df <- data.frame(Gene = character(), ALT = character())
}

## Report indels in gene regions
write_tsv(
  indels.Feat.df,
  file.path(
    varDir,
    paste0(
      argv$samplegroup,
      "genefiltIndels.tsv"
    )
  )
)

#### Get gene details for SNPs ####
snpCDSAttr <- pf_featuresNovar[
  findOverlaps(GRanges(AApred), GRanges(pf_featuresNovar),
    select = "last"
  )
]

snps.Feat.df <- data.frame(
  seqnames = seqnames(AApred),
  pos = start(AApred),
  ID = getGffAttribute(snpCDSAttr, "ID") |> unlist()
)
baseIDEvents <- str_replace(
  snps.Feat.df$ID, "^.*P", "P"
) |>
  str_remove("[\\.\\-].*$")
geneDetail <- map(baseIDEvents, function(ID) {
  genegff <- pf_featuresNovar[getGffAttribute(pf_featuresNovar, "ID") == ID]
  data.frame(
    GeneName = getGffAttribute(genegff, "Name"),
    Description = getGffAttribute(genegff, "description"),
    GeneID = ID
  )
}) |>
  list_rbind()
snps.Feat.df <- cbind(snps.Feat.df, geneDetail) |>
  mutate(
    Gene = case_when(
      is.na(GeneName) ~ Description |>
        str_replace_all(pattern = "\\+", replacement = " "),
      TRUE ~ GeneName
    ),
    GeneName = NULL, Description = NULL, ID = NULL
  )

SNPdetails <- left_join(
  SNPdetails, snps.Feat.df
)

write_tsv(
  SNPdetails,
  file.path(
    varDir, paste0(argv$samplegroup, "nonsynSNPs", ".tsv")
  )
)

## Write summary table
events.stats.df <- data.frame(
  initialnum = length(samplevcf),
  filtnum = length(eventFilt),
  majoritySom = length(majsom),
  SNP_CDS = length(snpCDS),
  SNP_nonsyn = length(AApred),
  Indel_gene = length(indelGene)
)
write_tsv(
  events.stats.df,
  file.path(
    varDir, paste0(
      argv$samplegroup, critSamplesSom,
      "plusstats.Qcrit", argv$QUALcrit, ".tsv"
    )
  )
)
