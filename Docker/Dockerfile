FROM continuumio/miniconda3:24.3.0-0

# install R4.4 from r repos since conda doesn't yet have 4.4
RUN apt-get update && apt-get install gpg -y
RUN echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" >> /etc/apt/sources.list && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys B8F25A8A73EACF41 && apt-get update && apt-get install r-base-dev -y && apt-get clean

# install R package external dependencies
RUN apt-get install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcurl4-openssl-dev libssl-dev -y && apt-get clean
# install R + bioconductor packages
RUN MAKE="make -j" Rscript -e 'install.packages(c("BiocManager", "devtools", "tidyverse", "argparser", "stringr", "readr", "scales", "gridExtra", "stringdist", "R6", "ggplot2"), repos="https://cloud.r-project.org/")'
RUN MAKE="make -j" Rscript -e 'BiocManager::install(c("biobase", "biocparallel", "QDAseq", "VariantAnnotation", "GenomeInfoDb", "SparseArray", "GenomicAlignments", "GenomeIntervals", "StructuralVariantAnnotation", "rtracklayer")); devtools::install_github("Bioconductor/BSgenome@e89b08b590039643f8fe9f030c1948233de4f50f")'

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
    multiqc=1.12.0 && conda clean -ay
