## Read stored copy numbers from copynumQDNAseq.R, and make plots
## of specified region of interest

suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())
pdf.options(useDingbats=FALSE)
suppressPackageStartupMessages(library(argparser))

#### Read command-line arguments ####
argp <- arg_parser(paste(
  "Read a scaled copy numbers data-frame file written by",
  "copynumQDNAseq and plot scaled CN for region of interest"
))
argp <- add_argument(argp, "--samplegroup",
                     help = "Group name of related samples. Required"
)
argp <- add_argument(argp, "--chrOI",
                     type = "character",
                     help = "Chromosome of region, as 2 digit integer",
                     default = "08"
)
argp <- add_argument(argp, "--startROI",
                     type = "numeric",
                     help = "start of ROI in kilobases",
                     default = 395
)
argp <- add_argument(argp, "--endROI",
                     type = "numeric",
                     help = "end of ROI in kilobases",
                     default = 435
)
argp <- add_argument(argp, "--parentId",
                     help = "parent name of related samples. Required"
)
argp <- add_argument(argp, "--bin_in_kbases",
                     default = "1",
                     help = "bin size in kbp"
)

argv <- parse_args(argp)

bin_in_kbases <- argv$bin_in_kbases
## Parent samples and samples of interest
parentId <- argv$parentId 
groupId <- argv$samplegroup

#### Read stored results ####
## Only CN scaled relative to parent used for zoomed-in plots
scaled_df <- tryCatch(
  readRDS( paste0(groupId, ".CN_compare_df_",
                  bin_in_kbases, "k", ".rds")
  ),
  error = function(e){
    stop(paste(e, "when reading copy number file",
               paste0(groupId, ".CN_compare_df_",
                      bin_in_kbases, "k", ".rds"))
    )
  }
)

## Copynum files are made with single parent,
## with sample name matching parentId in groupkey
sampleL <- colnames( scaled_df )[which( 
    str_detect( colnames( scaled_df ), parentId ) ) ] %>%
    str_remove( paste( ' vs.', parentId ) )

#### Define functions ####
samplecolrs <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(sampleL))
plotZoomedROI <- function(bin_df, startkb = 395, endkb = 435, chro = "08") {
    bin_df$chrom <- sapply(bin_df$chrom, function(x){
        strsplit(x, split='_')[[1]][2] } ) # shorten chromosome names
    chromOI <- pivot_longer(bin_df[bin_df$chrom==chro, ],
                            cols = !c("chrom", "range", "start", "end"), 
                            names_to = "sample", values_to = "copynum")
    chromOI$sample <-  sub(paste(" vs.", parentId), "", chromOI$sample)
    # set undefined copynum (presumed 0/0) to 1, and mark with a different shape
    chromOI$naCN <- is.na(chromOI$copynum)
    chromOI$copynum[is.na(chromOI$copynum)] <- 1  
    roi <- chromOI[which(chromOI$start > startkb*1000 & 
                             chromOI$end < endkb*1000), ]
    ymax <- ceiling(max(roi$copynum)) 
    ggplot(roi, aes(x = start/1000, xend=end/1000, 
                    y = copynum, yend=copynum) ) + 
      geom_point( aes(color = sample, shape = naCN) ) + 
      scale_shape_manual(values = c(19, 17)) +
      facet_grid(sample ~ ., scales = "free_y") + 
      theme(legend.position = "none") + 
      scale_colour_manual(values=samplecolrs) + 
      scale_y_continuous(limits=c(0, ymax),
                         breaks=seq(0, ymax, 2)) + 
      labs( x= paste( "Position in chromosome", chro, "(kb)" ),
            y="Scaled Copy number")
}

#### Zoom to regions of interest ####
panelplot <- plotZoomedROI(
  scaled_df,
  chro = ifelse( argv$chrOI %in% seq(1,9) |> as.character(),
                 ## check for single digit in chrom name
                 paste0("0", argv$chrOI),
                 argv$chrOI),
  startkb = argv$startROI,
  endkb = argv$endROI
) +
  labs(title = paste(
    bin_in_kbases, "kb copy numbers zoomed in on region of interest"
  ))

ggsave(
  filename = file.path(
     paste0(groupId, ".", "CNcolourpanels_chr", argv$chrOI, "_roi.pdf")
  ),
  units = "mm", width = 180, height = 240
)
