GenomicBreaks
=============

Purpose
-------

_GenomicBreaks_ is a `R` package using _[Bioconductor](https://bioconductor.org/)_
libraries to analyse pairwise alignments of whole genomes in which the gene
order has been scrambled by evolution, like in the picture below that represents
the comparison of homologous chromosomes in two distantly related molds,
_N. crassa_ (chrIII) and _P. comata_ (chr7).

<center><img src="man/figures/plotApairOfChrs_Neu-2.png" alt="Comparison between Neurospora crassa chrIII / Podospora comata chr7 (rev-complemented)" width = "40%"/></center>

This package is especially designed to parse and process the alignment files
produced by the our [pairwise genome alignment
pipeline](https://github.com/oist/plessy_pairwiseGenomeComparison), but should
be capable to import output of other pipelines as well.

Installation
------------

### Install the package.

The following should work:

```
Rscript -e 'remotes::install_github("oist/GenomicBreaks", repos=BiocManager::repositories())'
```

Add `dependencies=TRUE` if you would like to install the packages needed to build the vignettes.

### How to create a Singularity container with GenomicBreaks installed.

See the [Singularity reciepe file](./Singularity.def)

```
sudo singularity build GenomicBreaks.sif Singularity.def
```

GenomicBreaks in brief:
-----------------------

A pairwise alignment of two genomes is loaded in `GBreaks` objects wrapping
the `GRanges` class.  Here is an example:

```
GBreaks object with 505 ranges and 2 metadata columns:
        seqnames      ranges strand |     score                   query
           <Rle>   <IRanges>  <Rle> | <numeric>               <GRanges>
    [1]     chrI  5860-10010      + |     10609  NC_047487.1:8723-12716
    [2]     chrI 11157-11782      + |      1253 NC_047487.1:13035-13680
    [3]     chrI 25371-26528      + |      3651 NC_047487.1:15380-16537
    [4]     chrI 26849-29699      + |      6130 NC_047487.1:17263-20175
    [5]     chrI 29937-30607      + |       767 NC_047487.1:20203-20866
    ...      ...         ...    ... .       ...                     ...
  [501]     chrM 67039-67870      + |      1441   NC_018044.1:7758-8439
  [502]     chrM 68205-68580      + |       830   NC_018044.1:8783-9180
  [503]     chrM 69178-76168      + |     14528  NC_018044.1:9650-16261
  [504]     chrM 77356-80022      + |      8066 NC_018044.1:53206-55865
  [505]     chrM 80919-85779      + |      6712 NC_018044.1:57885-61592
  -------
  seqinfo: 17 sequences (1 circular) from sacCer3 genome
```

See “_Get started_” on <https://oist.github.io/GenomicBreaks> for further details.
