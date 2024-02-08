## Read stored copy numbers from copynumQDNAseqParents.R, and make plots
## of nuclear genome

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library("Biobase"))
theme_set(theme_bw())
pdf.options(useDingbats=FALSE)

suppressPackageStartupMessages(library(gridExtra))  # to arrange the 14 chromosomes in 2 rows
suppressPackageStartupMessages(library(scales) )   # Need the scales package for log2 transform
suppressPackageStartupMessages(library(argparser))

#### Read command-line arguments ####
argp <- arg_parser(paste(
  "Read stored copy numbers, and make plots."
))
argp <- add_argument(argp, "--samplegroup",
                     help = "Group name of related samples. Required"
)
argp <- add_argument(argp, "--parentId",
                     help = "parent name of related samples. Required"
)
argp <- add_argument(argp, "--refDir",
                     help = "reference directory "
)
argp <- add_argument(argp, "--bin_in_kbases",
                     default = "5",
                     help = "bin size in kbp"
)
argp <- add_argument(argp, "--lowerbound_plot",
                     type = "numeric",
                     help = "minCN for whole-genome plot. Optional"
)
argp <- add_argument(argp, "--upperbound_plot",
                     type = "numeric",
                     help = "maxCN for whole-genome plot. Optional"
)

argv <- parse_args(argp)
refDir <- argv$refDir

bin_in_kbases <- argv$bin_in_kbases
## Parent samples and samples of interest
parentId <- argv$parentId 
groupId <- argv$samplegroup

#### Read stored results ####
copyNumsNormed <- tryCatch(
  readRDS( paste0( groupId, '.CN_df_', bin_in_kbases, 'k.rds')
                ),
  error <- function(e){
    stop(paste(e, "when reading copy number file",
               paste0( groupId, '.CN_df_', bin_in_kbases, 'k.rds') )
         )
  }
)
scaled_df <- tryCatch(
  readRDS( paste0(groupId, ".CN_compare_df_",
                  argv$bin_in_kbases, "k", ".rds")
  ),
  error <- function(e){
    stop(paste(e, "when reading copy number file",
               paste0(groupId, ".CN_compare_df_",
                      argv$bin_in_kbases, "k", ".rds"))
    )
  }
)
scaled_nuc <- dplyr::filter(scaled_df, !(chrom %in% nonNuc) ) 

## Copynum files are made with single parent,
## with sample name matching parentId in groupkey
sampleL <- colnames( scaled_nuc )[which( 
    str_detect( colnames( scaled_nuc ), parentId ) ) ] %>%
    str_remove( paste( ' vs.', parentId ) )

#### Define functions ####

plotLogRow <- function(binCounts, orgCol, bin_size, maxCN, minCN) {
    ## Plot chromosomes in single row, log-y scale
    ## ymax and ymin can be passed in or calculated from data
    binCounts$score <- binCounts[, orgCol]
    binCounts$shortChrom <- sapply(binCounts$chrom, function(x){
        strsplit(x, split='_')[[1]][2] } )
    binCounts$pos <- paste(binCounts$shortChrom, binCounts$start, sep=":") %>%
        factor(., levels = paste(binCounts$shortChrom, binCounts$start, sep=":") )
    chromlist <- unique(binCounts$shortChrom)
    if (missing(maxCN) | is.na(maxCN)) maxCN <- max(na.omit(binCounts$score))
    if (missing(minCN) | is.na(minCN)) minCN <- min(
      binCounts[which(binCounts$score > 0), 'score'])
    binCounts$score <- pmax(binCounts$score, minCN)
    binCounts$offscale <- binCounts$score > maxCN
    binCounts$score <- pmin(binCounts$score, maxCN)
    ggplot(binCounts , aes(pos, score, shape=offscale)) +
        geom_point(size=0.5) +
        scale_x_discrete(breaks=paste(chromlist, "1", sep=":"),
                         labels=paste('  ', chromlist) ) +
        labs(title = paste(orgCol, bin_size, "k bins Copy Numbers"),
             x='', y = '') +
        scale_y_continuous(trans=log2_trans(), limits=c(minCN, 1.1*maxCN),
                           breaks = trans_breaks("log2", function(x) 2^x),
                           labels = trans_format("log2", math_format(2^.x))) +
        scale_shape_manual(values = c(16, 18), guide = "none")
}

#### Make a plot for parent plus un-scaled samples in a single column ####
## Option to fix values for ymax and ymin to give better comparisons
nucCN <- copyNumsNormed %>%
    dplyr::filter( !(chrom %in% nonNuc) )

plotset <- arrangeGrob(
    grobs = lapply(c( parentId, sampleL[order(sampleL)] ), function(samplen)
    { plotLogRow(
      allCN, samplen, argv$bin_in_kbases,
      maxCN = argv$upperbound_plot,
      minCN = argv$lowerbound_plot
    ) + theme(title = element_text(size = 8),
              plot.margin = margin(0,4,0,4),
              axis.text = element_text(size = 6))
        }),
    ncol = 1)
plot(plotset)
ggsave( plotset, filename = file.path(paste0(groupId, "_copynums_", 
                    bin_in_kbases, "k.pdf") )
            , units = "mm", width = 297, height = 210 
        )
