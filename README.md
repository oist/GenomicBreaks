GenomicBreaks
=============

Purpose
-------

We are developing the `GenomicBreaks` package to parse and process the
alignment files produced by the our
[pariwise genome alignment pipeline](https://github.com/oist/plessy_pairwiseGenomeComparison).

Installation
------------

### Install the package.

The following should work:

```
Rscript -e 'remotes::install_github("oist/GenomicBreaks", repos=BiocManager::repositories())'
```

### How to install R 4.1, Rstudio and Bioconductor.

On a Debian/Ubuntu system, try this:

```
sudo apt install r-base
sudo apt install pandoc qpdf texlive # For vignette builds and package checks
sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev libfftw3-dev libtiff-dev libgsl-dev
sudp atp install libfontconfig1-dev libharfbuzz-dev libfribidi-dev # For pkgdown
sudo apt install git bash-completion
sudo apt install libgl1 libnss3 libasound2 libxdamage1
wget https://s3.amazonaws.com/rstudio-ide-build/desktop/bionic/amd64/rstudio-1.4.1717-amd64.deb
sudo apt --fix-broken install ./rstudio-1.4.1717-amd64.deb
Rscript -e 'install.packages("BiocManager")'
Rscript -e 'BiocManager::install("GenomicRanges")'
Rscript -e 'install.packages("tidyverse")'
Rscript -e 'install.packages("devtools")' 
Rscript -e 'install.packages("remotes")'
```

GenomicBreaks in brief:
-----------------------

A pairwise alignment of two genomes is loaded in `GBreaks` objects wrapping
the `GRanges` class.  Here is an example:

```
GRanges object with 11 ranges and 2 metadata columns:
       seqnames            ranges strand |                  query     score
          <Rle>         <IRanges>  <Rle> |              <GRanges> <integer>
   [1]     chr1 11234702-11242676      + |   Chr1:7691083-7699807      7975
   [2]     chr1 11242764-11256464      - |   Chr1:7714830-7733148     13701
   [3]     chr1 11256821-11271214      - |   Chr1:7699877-7713142     14394
   [4]     chr1 11271261-11272159      - |   Chr1:7975442-7976321       899
   [5]     chr1 11272246-11274272      + |   Chr1:7686802-7688942      2027
   [6]     chr1 11275227-11276200      - |   Chr1:7491169-7492136       974
   [7]     chr1 11276902-11281111      - |   Chr1:7850371-7855204      4210
   [8]     chr1 11281154-11281731      + |  PAR:10891068-10891635       578
   [9]     chr1 11281946-11288799      + |   Chr2:9359434-9367027      6854
  [10]     chr1 11288839-11299743      - | Chr1:10912857-10921537     10905
  [11]     chr1 11300902-11301564      - |   Chr1:9597979-9599493       663
  -------
  seqinfo: 19 sequences from OKI2018_I69 genome
```

See <https://oist.github.io/GenomicBreaks> for further details.
