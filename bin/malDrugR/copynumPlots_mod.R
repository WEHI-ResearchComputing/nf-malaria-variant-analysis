## Read stored copy numbers from copynumQDNAseqParents.R, and make plots
## Basically same as version in Madeline 2022

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
                     default = "1",
                     help = "bin size in kbp"
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

plotLogRow <- function(binCounts, orgCol, bin_size, maxCN, minCN, newtitle) {
    ## Plot chromosomes in single row, log-y scale
    ## ymax and ymin can be passed in or calculated from data
    binCounts$score <- binCounts[, orgCol]
    binCounts$shortChrom <- sapply(binCounts$chrom, function(x){
        strsplit(x, split='_')[[1]][2] } )
    binCounts$pos <- paste(binCounts$shortChrom, binCounts$start, sep=":") %>%
        factor(., levels = paste(binCounts$shortChrom, binCounts$start, sep=":") )
    chromlist <- unique(binCounts$shortChrom)
    if (missing(maxCN)) maxCN <- max(na.omit(binCounts$score))
    if (missing(minCN)) minCN <- min(
        binCounts[which(binCounts$score > 0), 'score'])
    if (missing(newtitle)) newtitle <- 
        paste( orgCol, bin_size, "k bins Copy Numbers" )
    binCounts$score <- pmax(binCounts$score, minCN)
    binCounts$offscale <- binCounts$score > maxCN
    binCounts$score <- pmin(binCounts$score, maxCN)
    ggplot(binCounts , aes(pos, score, shape=offscale)) +
        geom_point(size=0.5) +
        scale_x_discrete(breaks=paste(chromlist, "1", sep=":"),
                         labels=paste('  ', chromlist) ) +
        labs(title = newtitle, x='', y = '') +
        scale_y_continuous(trans=log2_trans(), limits=c(minCN, 1.1*maxCN),
                           breaks = trans_breaks("log2", function(x) 2^x),
                           labels = trans_format("log2", math_format(2^.x))) +
        scale_shape_manual(values = c(16, 18), guide = "none")
}
deemph1 <- function(cn){
    ifelse( is.na( cn ), 0,
            abs( plogis( cn, location = 1, scale = 0.5 ) -0.5 )*2
    )
    }
plotWholeCN <- function(bin_df, bin_size, maxCN){
    # plot copynumbers in a single panel, coloured by Sample
    # cap the number of bins at 5000, plus start (for labels)
    bin_df$chrom <- sapply(bin_df$chrom, function(x){
        strsplit(x, split='_')[[1]][2] } ) # shorten chromosome names
    bin_df$pos <- paste(bin_df$chrom, bin_df$start, sep=":") %>%
        factor(., levels = paste(bin_df$chrom, bin_df$start, sep=":") )
    bin_df <- reshape2::melt(
        bin_df, id.vars = c("chrom", "range", "start", "end", "pos"),
        variable.name = "sample", value.name = "copynum") 
    bin_df$sample <-  sub(paste( " vs.", parentId), "",bin_df$sample)
    chromlist <- unique(bin_df$chrom)
    if (missing(maxCN)) maxCN <- ceiling( max(na.omit(bin_df$copynum)) )
    smaller_df <- 
        rbind( filter(bin_df,  start == 1),
               slice_sample( bin_df, n=5000, 
                             weight_by = deemph1(bin_df$copynum) + 0.1)
        )
    # x-axis labelling requires all chromosomes have a row with start==1 
    # Will fail if na.omit used
    ggplot(smaller_df, aes(x=pos, y=copynum, colour=sample)) +
        geom_point(size=0.5
                   ) +
        scale_colour_brewer(palette="Set1") +
        scale_x_discrete( drop = FALSE,
                          breaks=paste(chromlist, "1", sep=":"),
                         labels=paste('  ', chromlist) ) +
        scale_y_continuous(limits=c(0, maxCN),
                           breaks=seq(0, maxCN, 2)) +
        labs(title= paste(bin_size,"kb copy numbers"), x= "Chromosome", y="Relative copy numbers") 
}

plotZoomedROI <- function(bin_df, startkb = 395, endkb = 435, chro = "08") {
    bin_df$chrom <- sapply(bin_df$chrom, function(x){
        strsplit(x, split='_')[[1]][2] } ) # shorten chromosome names
    chromOI <- pivot_longer(bin_df[bin_df$chrom==chro, ],
                            cols = !c("chrom", "range", "start", "end"), 
                            names_to = "sample", values_to = "copynum")
    chromOI$sample <-  sub(paste(" vs.", parentId), "",chromOI$sample)
    # set undefined copynum (presumed 0/0) to 1, and mark with a different line-type
    chromOI$naCN <- is.na(chromOI$copynum)
    chromOI$copynum[is.na(chromOI$copynum)] <- 1  
    roi <- chromOI[which(chromOI$start > startkb*1000 & 
                             chromOI$end < endkb*1000), ]
    ymax <- ceiling(max(roi$copynum)) 
    ggplot(roi, aes(x = start/1000, xend=end/1000, 
                    y = copynum, yend=copynum) ) + 
        geom_segment( aes(color = sample, linetype = naCN) ) + 
        facet_grid(sample ~ ., scales = "free_y") + 
        theme(legend.position = "none") + 
        scale_colour_brewer(palette="Dark2") + 
        scale_y_continuous(limits=c(0, ymax),
                           breaks=seq(0, ymax, 2)) + 
        labs( x= paste( "Position in chromosome", chro, "(kb)" ), y="Scaled Copy number")
}

## Make a plot for parent plus un-scaled samples in a single column
## Option to fix values for ymax and ymin to give better comparisons
nucCN <- copyNumsNormed %>%
    dplyr::filter( !(chrom %in% nonNuc) )

plotset <- arrangeGrob(
    grobs = lapply(c( parentId, sampleL[order(sampleL)] ), function(samplen)
    { plotLogRow(nucCN, paste(samplen) , argv$bin_in_kbases
                 , maxCN = 2^5, minCN = 2^-2
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

#### Zoom to regions of interest ####
## PF samples initial ROI is 08: 410-435; other is 14: 1182-1219
chrname <- "08"
roiStartkb <- 300; roiEndkb <- 600 
panelplot <- plotZoomedROI(
    scaled_df, startkb = roiStartkb 
    , endkb = roiEndkb 
    , chro = chrname )
panelplot + labs(title=paste(
    bin_in_kbases, "kb copy numbers zoomed in on region of interest"))

ggsave(filename = file.path(
            paste0(groupId, ".", "CNcolourpanels_chr", chrname, "_roi.pdf") )
            , units = "mm", width = 180, height = 240
        )
