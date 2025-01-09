FROM continuumio/miniconda3:24.3.0-0

# install remaining packages with conda
RUN conda update --all -c conda-forge --override-channels && \
    conda install -c bioconda -c conda-forge --override-channels -y \
    bwa=0.7.17 \
     && conda clean -ay
