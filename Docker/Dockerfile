FROM ubuntu:22.04

# get conda
RUN apt-get update && apt-get install wget -y
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh && bash /opt/miniconda.sh -b -u -p /opt/miniconda3 && rm /opt/miniconda.sh
ENV PATH=/opt/miniconda3/bin:$PATH
# install packages
RUN conda install -c bioconda -c conda-forge --override-channels -y \
    bcftools=1.13 \
    bwa=0.7.17 \
    samtools=1.16.1 \
    picard=2.26.11 \
    fastqc=0.11.8-0 \
    multiqc=1.12.0 \
    r-base=4.3.2 \
    r-devtools \
    r-tidyverse \
    r-argparser \
    r-stringr \
    r-readr \
    r-scales \
    r-gridextra \
    r-stringdist \
    r-r6 \
    r-ggplot2 \
    bioconda::bioconductor-biobase \
    bioconda::bioconductor-biocparallel \
    bioconda::bioconductor-bsgenome \
    bioconda::bioconductor-qdnaseq \
    bioconda::bioconductor-variantannotation \
    bioconda::bioconductor-genomeinfodb \
    bioconda::bioconductor-sparsearray \
    bioconda::bioconductor-genomicalignments \
    bioconda::bioconductor-genomeintervals \
    bioconda::bioconductor-structuralvariantannotation \
    bioconda::bioconductor-rtracklayer \
    gridss=2.13.2 && conda clean -ay

COPY entrypoint.sh /usr/local/bin/entrypoints.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]