# Get genomic sequences

Extract the sequences of the genomic ranges from a
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
or a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object.

## Usage

``` r
getSeq(x, ...)

# S4 method for class 'GRanges'
getSeq(x, ...)

# S4 method for class 'GBreaks'
getSeq(x, ...)
```

## Arguments

- x:

  A `GBreaks` or a `GRanges` object.

- ...:

  Extra arguments (not used, but present for compatibility).

## Value

A
[`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object containing the extracted sequence(s).

## See also

See also the
[`Biostrings::getSeq`](https://rdrr.io/pkg/Biostrings/man/getSeq.html)
function.

Other Bioconductor API functions:
[`GBreaks-class`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md),
[`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md),
[`range_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/range.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`subsetByOverlaps_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/subsetByOverlaps.md)

## Author

Charles Plessy

## Examples

``` r
Scerevisiae <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
getSeq(GRanges("chrI:1-20", seqinfo = seqinfo(Scerevisiae)))
#> DNAStringSet object of length 1:
#>     width seq
#> [1]    20 CCACACCACACCCACACACC
```
