FROM rocker/r-ver:4.4.3

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    cmake \
    git \
    wget \
    ca-certificates \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libblas-dev \
    liblapack-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Conda via Miniforge (plus propre pour conda-forge)
ENV CONDA_DIR=/opt/conda
ENV PATH=/opt/conda/bin:$PATH
ENV R_HOME=

RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm -f /tmp/miniforge.sh

WORKDIR /app
COPY requirements_py.txt /app/requirements_py.txt

# Créer l'env conda
RUN conda create -y -n BulkTools -c conda-forge \
    r-base=4.4.3 \
    python=3.12 \
    pip

# Paquets R fondamentaux installés via conda-forge
RUN conda run -n BulkTools conda install -y -c conda-forge \
    r-shiny \
    r-shinyfiles \
    r-fs \
    r-ggrepel \
    r-glmnet \
    r-survminer \
    r-survival \
    r-optparse \
    r-ggplot2 \
    r-dplyr \
    r-dt \
    r-bslib \
    r-jsonlite \
    r-data.table \
    r-matrix

# Packages Python
RUN conda run -n BulkTools pip install --no-cache-dir -r /app/requirements_py.txt

# CRAN packages
RUN conda run -n BulkTools R -q -e "install.packages(c('tidyestimate'), repos='https://cloud.r-project.org', Ncpus=max(1, parallel::detectCores()-1))"
    
# Download R libraries through conda bioconductor
RUN conda install -y -n BulkTools \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority \
    r-tidyverse \
    r-reticulate \
    bioconductor-tximport \
    bioconductor-genomeinfodb \
    bioconductor-annotationdbi \
    bioconductor-fgsea \
    bioconductor-gsva \
    bioconductor-deseq2 \
    bioconductor-limma \
    bioconductor-outrider \
    bioconductor-txdb.hsapiens.ucsc.hg19.knowngene \ 
    bioconductor-org.hs.eg.db


# test
RUN conda run -n BulkTools R -q -e "library(shiny); library(shinyFiles); library(fs); library(optparse); library(tidyverse); cat('Docker image 1/2 OK\\n')"
RUN conda run -n BulkTools R -q -e "library(org.Hs.eg.db); library(tidyestimate); library(tximport); library(GenomeInfoDb); library(AnnotationDbi); cat('Docker image 2/3 OK\\n')"
RUN conda run -n BulkTools R -q -e "library(fgsea); library(GSVA); library(DESeq2); library(jsonlite); library(data.table); library(Matrix); library(OUTRIDER); cat('Docker image 3/3 OK\\n')"

COPY demarreur_app.R /app/demarreur_app.R
COPY app/ /app/
COPY REF_DATA/ /app/REF_DATA/


EXPOSE 5288

ENTRYPOINT ["/opt/conda/envs/BulkTools/bin/Rscript", "/app/demarreur_app.R"]