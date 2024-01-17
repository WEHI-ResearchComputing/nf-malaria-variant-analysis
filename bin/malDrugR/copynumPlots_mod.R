## Read stored copy numbers from copynumQDNAseqParents.R, and make plots
## Basically same as version in Madeline 2022

library(tidyverse)
library("QDNAseq"); library("Biobase")
theme_set(theme_bw())
pdf.options(useDingbats=FALSE)
## This should remove need for 'useDingbats=FALSE' in ggsave commands.
# Alternative solution is to edit fonts used by Illustrator

library(gridExtra)  # to arrange the 14 chromosomes in 2 rows
library(scales)     # Need the scales package for log2 transform

PAPDIR <- "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects" 

pfCurrRefs <- data.frame(
    strain = c( '3D7', 'Dd2' ),
    version = c( '52', '57' )
)
ref <- pfCurrRefs[1,]
refDir <- file.path( 
    PAPDIR, "reference_genomes", "plasmodium",
    paste0( "PlasmoDB-", ref$version, "_Pfalciparum", ref$strain) )
WORKDIR <- "/vast/scratch/users/penington.j/malaria"
alignDir <- file.path( WORKDIR, "alignbwa" )
cnDir <- file.path( WORKDIR, "variants", "copynumfiles")
setwd( WORKDIR )
plotDir <- file.path( WORKDIR, "variants", "copynumPlots")
if( !dir.exists( plotDir ) ){ dir.create( plotDir ) }

bin_in_kbases <- "5"
## Parent samples and samples of interest
parentName <- 'S107' 
analysis <- 'S108vsS107'
strain <- # read sample_key

#### Read stored results ####
copyNumsNormed <- readRDS( file.path(
    cnDir, paste0( analysis, '_cn_', bin_in_kbases, 'k.rds')
) )
scaled_df <- readRDS(
    file.path(cnDir, 
              paste0(analysis, ".CN_compare_df_", 
                     bin_in_kbases, "k", ".rds"))) 
## Remove Apical and Mitochondrial genome bins from data if present, 
## and re-read data later for mito inspection. Genomes too small to show on 
## same plot as nuclear
if(ref$strain == 'Dd2') {nonNuc <- c( 'PfDd2_API', 'PfDd2_MT')
} else {nonNuc <- c('Pf3D7_API_v3', 'Pf3D7_MIT_v3') }

scaled_df <- dplyr::filter(scaled_df, !(chrom %in% nonNuc) ) 

## Don't Read sample names from data frame any more. Use sample_key. 
## No strain name, just analysis name, or optional input
sampleL <- colnames( scaled_df )[which( 
    str_detect( colnames( scaled_df ), strain ) ) ] %>%
    str_remove( paste( ' vs.', parentName ) )
## - str_detect(..., strain) doesn't work when not all have strain as prefix

#### Define functions ####
convertQDNAtoDF <- function(qobject) {
    if ("counts" %in% assayDataElementNames(qobject) ) {
        scoreDF <- as.data.frame( assayDataElement(qobject, "counts"))
    } else if ("copynumber" %in% assayDataElementNames(qobject) ) {
        scoreDF <- as.data.frame( assayDataElement(qobject, "copynumber"))
    }
    scoreDF$chrom <- sapply(strsplit(rownames(scoreDF), ':'), function(x) (unlist(x)[1]) )
    scoreDF$range <- sapply(strsplit(rownames(scoreDF), ':'), function(x) (unlist(x)[2]) )
    scoreDF$start <- sapply(strsplit(scoreDF$range, '-'), function(x) (as.numeric(unlist(x)[1]) ) )
    scoreDF$end <- sapply(strsplit(scoreDF$range, '-'), function(x) (as.numeric(unlist(x)[2]) ) )
    
    return(scoreDF)
}

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
    bin_df$sample <-  sub(paste( " vs.", parentName), "",bin_df$sample)
    chromlist <- unique(bin_df$chrom)
    if (missing(maxCN)) maxCN <- ceiling( max(na.omit(bin_df$copynum)) )
    smaller_df <- 
        rbind( filter(bin_df,  start == 1 ),
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
    chromOI$sample <-  sub(paste(" vs.", parentName), "",chromOI$sample)
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
allCN <- convertQDNAtoDF(copyNumsNormed) %>%
    dplyr::filter( !(chrom %in% nonNuc) )

plotset <- arrangeGrob(
    grobs = lapply(c( parentName, sampleL[order(sampleL)] ), function(samplen)
    { plotLogRow(allCN, paste(samplen ) , bin_in_kbases
                 , maxCN = 2^5, minCN = 2^-2
    ) + theme(title = element_text(size = 8),
              plot.margin = margin(0,4,0,4),
              axis.text = element_text(size = 6))
        }),
    ncol = 1)
plot(plotset)
ggsave( plotset, filename = file.path(
    plotDir, paste0(analysis, "_copynums_", 
                    bin_in_kbases, "k.pdf") )
    , units = "mm", width = 297, height = 210 )

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
    plotDir, paste0(strain, ".", "CNcolourpanels_chr", chrname, "_roi.pdf") )
    , units = "mm", width = 180, height = 240
    # , units = "mm", width = 70, height = 70
)
