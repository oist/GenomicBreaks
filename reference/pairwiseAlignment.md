# Pairwise alignment of genomic ranges

Retrieves DNA sequence from
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
or
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
objects that are properly related to a `BSgenome` package, and aligns
them with the `pwalign::pairwiseAlignment()` function.

## Usage

``` r
pairwiseAlignment(pattern, subject, ...)

# S4 method for class 'GRanges,GRanges'
pairwiseAlignment(pattern, subject, ...)

# S4 method for class 'GBreaks,ANY'
pairwiseAlignment(pattern, subject, ...)
```

## Arguments

- pattern:

  A `GBreaks` or a `GRanges` object

- subject:

  A `GBreaks` or a `GRanges` object

- ...:

  Additional arguments passed to `pairwiseAlignment`.

## Value

Returns a
[`pwalign::PairwiseAlignments`](https://rdrr.io/pkg/pwalign/man/PairwiseAlignments-class.html)
object.

## See also

Other Bioconductor API functions:
[`GBreaks-class`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md),
[`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md),
[`range_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/range.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`subsetByOverlaps_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/subsetByOverlaps.md)

## Author

Charles Plessy

## Examples

``` r
Scerevisiae <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
# Very arbitrary example
gb <- GBreaks( target = GRanges("chrI: 1-20", seqinfo = seqinfo(Scerevisiae))
             , query  = GRanges("chrI:21-40", seqinfo = seqinfo(Scerevisiae)))
pairwiseAlignment(gb)
#> Error in .call_fun_in_pwalign("pairwiseAlignment", ...): pairwiseAlignment() has moved from Biostrings to the pwalign package, and is
#>   formally defunct in Biostrings >= 2.77.1. Please call
#>   pwalign::pairwiseAlignment() to get rid of this error.
```
