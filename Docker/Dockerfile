FROM mambaorg/micromamba:git-79bf814-ubuntu24.04

# install R, R package external dependencies and other tools with conda
# environment-linux-64.lock is generated from environment.yml, regenerate it after editing that file with:
#   conda-lock lock -f environment.yml -p linux-64 --kind explicit --filename-template 'environment-{platform}.lock'
COPY --chown=$MAMBA_USER:$MAMBA_USER environment-linux-64.lock /tmp/environment.lock
RUN micromamba install -y -n base -f /tmp/environment.lock && \
    micromamba clean -ay

# activate the base environment for the RUN commands below, so they can find Rscript, wget etc
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# install R + bioconductor packages
RUN MAKE="make -j" Rscript -e 'install.packages(c("BiocManager", "devtools", "tidyverse", "argparser", "stringr", "readr", "scales", "gridExtra", "stringdist", "R6", "ggplot2"), repos="https://cloud.r-project.org/")'
RUN MAKE="make -j" Rscript -e 'BiocManager::install(c("biobase", "biocparallel", "QDNAseq", "VariantAnnotation", "GenomeInfoDb", "SparseArray", "GenomicAlignments", "genomeIntervals", "StructuralVariantAnnotation", "rtracklayer")); devtools::install_github("Bioconductor/BSgenome@e89b08b590039643f8fe9f030c1948233de4f50f", upgrade = FALSE)'

COPY BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz /tmp/BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz
COPY BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz /tmp/BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz
COPY BSgenome.PfalciparumNF54iGP.tar.gz /tmp/BSgenome.PfalciparumNF54iGP.tar.gz
RUN R CMD INSTALL /tmp/BSgenome.Pfalciparum3D7.PlasmoDB.52.tar.gz &&  \
    R CMD INSTALL /tmp/BSgenome.PfalciparumDd2.PlasmoDB.57.tar.gz &&  \
    R CMD INSTALL /tmp/BSgenome.PfalciparumNF54iGP.tar.gz

# install picard and gridss manually since conda versions depend on R
# root is needed to write to /opt and /usr/local/bin
USER root
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar -O /opt/picard.jar && echo '#!/bin/bash' > /usr/local/bin/picard && echo 'java -jar /opt/picard.jar "$@"' >> /usr/local/bin/picard && chmod +x /usr/local/bin/picard
WORKDIR /tmp
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz && mkdir /opt/gridss && cd /opt/gridss && tar xzvf /tmp/gridss-2.13.2.tar.gz && rm /tmp/gridss-2.13.2.tar.gz
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar -O /opt/gridss/gridss.jar
USER $MAMBA_USER
ENV PATH="/opt/gridss:$PATH"
ENV GRIDSS_JAR=/opt/gridss/gridss.jar
