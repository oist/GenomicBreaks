# The GenomicBreaks class

The `GBreaks` class is a simple wrapper to the
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
class.

## Details

Aligned sequences of the *target* genome are represented as the main
ranges of the `GRanges` object, and their counterparts in the *query*
genome are represented as a `GRanges` object sored in the the metadata
column `query`.

## See also

Other Bioconductor API functions:
[`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md),
[`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md),
[`range_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/range.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`subsetByOverlaps_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/subsetByOverlaps.md)

## Examples

``` r
gb <- load_genomic_breaks(system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"))
gb
#> GBreaks object with 505 ranges and 2 metadata columns:
#>         seqnames        ranges strand |     score                       query
#>            <Rle>     <IRanges>  <Rle> | <numeric>                   <GRanges>
#>     [1]     chrM       191-813      + |       984     NC_018044.1:61854-62304
#>     [2]     chrM     6313-9737      - |      9007     NC_018044.1:62854-65728
#>     [3]     chrM   12828-13986      + |      2733     NC_018044.1:22453-23578
#>     [4]     chrM   16435-18988      + |     12413     NC_018044.1:23579-26133
#>     [5]     chrM   20505-20650      + |       774     NC_018044.1:26134-26279
#>     ...      ...           ...    ... .       ...                         ...
#>   [501]  chrXIII 769669-808526      + |    152331   NC_047499.1:762339-801278
#>   [502]  chrXIII 809285-837102      + |    115538   NC_047499.1:801658-829480
#>   [503]  chrXIII 837914-912582      + |    303512   NC_047499.1:830467-905221
#>   [504]  chrXIII 912605-917242      + |     13091   NC_047499.1:905399-910082
#>   [505]  chrXIII 919705-922192      + |      6580 NC_047501.1:1054871-1057371
#>   -------
#>   seqinfo: 17 sequences from an unspecified genome; no seqlengths
gb$query
#> GRanges object with 505 ranges and 0 metadata columns:
#>            seqnames          ranges strand
#>               <Rle>       <IRanges>  <Rle>
#>     [1] NC_018044.1     61854-62304      *
#>     [2] NC_018044.1     62854-65728      *
#>     [3] NC_018044.1     22453-23578      *
#>     [4] NC_018044.1     23579-26133      *
#>     [5] NC_018044.1     26134-26279      *
#>     ...         ...             ...    ...
#>   [501] NC_047499.1   762339-801278      *
#>   [502] NC_047499.1   801658-829480      *
#>   [503] NC_047499.1   830467-905221      *
#>   [504] NC_047499.1   905399-910082      *
#>   [505] NC_047501.1 1054871-1057371      *
#>   -------
#>   seqinfo: 17 sequences from an unspecified genome; no seqlengths
```
