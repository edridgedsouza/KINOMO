#! /usr/bin/env bash
#  mamba is used for speed but conda also works here

mamba install r-base=4.3.1 r-cluster r-cowplot r-gplots r-nmf \
    r-pkgmaker r-rngtools r-seurat r-tidyverse \
    -c r -c conda-forge -c bioconda 
