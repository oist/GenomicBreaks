FROM bioconductor/bioconductor_docker:RELEASE_3_13
LABEL authors="charles.plessy@oist.jp" \
      description="GenomicBreaks package in the Bioconductor Docker image"

RUN R CMD INSTALL .
