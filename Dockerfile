FROM bioconductor/bioconductor_docker:RELEASE_3_13

LABEL authors="charles.plessy@oist.jp" \
      description="GenomicBreaks package in the Bioconductor Docker image"

RUN apt update && \
	apt install -y pandoc \
	qpdf \
	texlive \
	libxml2 \
	libcurl4-openssl-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	git \
	bash-completion \
	libgl1 \
	libnss3 \
	libasound2 \
	libxdamage1 \
	libbz2-dev \
	liblzma-dev \
	libfftw3-dev

RUN Rscript -e 'install.packages("BiocManager")' && \
	Rscript -e 'install.packages("tidyverse")' && \
	Rscript -e 'install.packages("devtools")' && \
	Rscript -e 'install.packages("remotes")' && \
	Rscript -e 'BiocManager::install("BSgenome")'

RUN Rscript -e 'remotes::install_github("oist/GenomicBreaks", repos=BiocManager::repositories(), dependencies=TRUE)' 
