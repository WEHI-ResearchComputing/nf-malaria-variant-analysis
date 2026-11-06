FROM continuumio/miniconda3:26.7.1

# install R package external dependencies
RUN apt-get update && apt-get install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcurl4-openssl-dev libssl-dev libwebp-dev pkg-config gnupg r-base-dev -y && apt-get clean
# install R + bioconductor packages
RUN MAKE="make -j" Rscript -e 'install.packages(c("BiocManager", "devtools", "tidyverse", "argparser", "stringr", "readr", "scales", "gridExtra", "stringdist", "R6", "ggplot2"), repos="https://cloud.r-project.org/")'
RUN MAKE="make -j" Rscript -e 'BiocManager::install(c("BSgenome", "biobase", "biocparallel", "QDNAseq", "VariantAnnotation", "GenomeInfoDb", "SparseArray", "GenomicAlignments", "genomeIntervals", "StructuralVariantAnnotation", "rtracklayer"))'

COPY BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz /tmp/BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz
COPY BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz /tmp/BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz
COPY BSgenome.PfalciparumNF54iGP.tar.gz /tmp/BSgenome.PfalciparumNF54iGP.tar.gz
RUN R CMD INSTALL /tmp/BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz &&  \
    R CMD INSTALL /tmp/BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz &&  \
    R CMD INSTALL /tmp/BSgenome.PfalciparumNF54iGP.tar.gz

# install picard and gridss manually since conda versions depend on R
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar -O /opt/picard.jar && echo '#!/bin/bash' > /usr/local/bin/picard && echo 'java -jar /opt/picard.jar "$@"' >> /usr/local/bin/picard && chmod +x /usr/local/bin/picard
WORKDIR /tmp
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz && mkdir /opt/gridss && cd /opt/gridss && tar xzvf /tmp/gridss-2.13.2.tar.gz && rm /tmp/gridss-2.13.2.tar.gz
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar -O /opt/gridss/gridss.jar
ENV PATH="/opt/gridss:$PATH"
ENV GRIDSS_JAR=/opt/gridss/gridss.jar


# install remaining packages with conda
RUN conda update --all -c conda-forge --override-channels && \
    conda install -c bioconda -c conda-forge --override-channels -y \
    bcftools=1.13 \
    bwa=0.7.17 \
    samtools=1.16.1 \
    fastqc=0.11.8-0 \
    mosdepth=0.3.8 \
    multiqc=1.29 \
    && conda clean -ay
