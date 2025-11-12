# Gaps between ranges

Utility function that runs
[`GenomicRanges::gaps()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)
and then cleans its output by removing strandless ranges as well as the
last range that represents the sequence between the end of the input
object and the end of its sequence levels.

## Usage

``` r
cleanGaps(gr)
```

## Arguments

- gr:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  or a
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

## Value

Returns a strandless
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object representing the gaps between the ranges of the input.

## Note

If you find replacement provided by a package that we already import,
please let me know in a GitHub issue or pull request.

## See also

Other Reducing functions:
[`get_bps()`](https://oist.github.io/GenomicBreaks/reference/get_bps.md),
[`wholeGenomeClassification()`](https://oist.github.io/GenomicBreaks/reference/wholeGenomeClassification.md)

## Author

Charles Plessy

## Examples

``` r
cleanGaps(exampleColinear)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA   151-250      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
GenomicRanges::gaps(exampleColinear)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA      1-99      +
#>   [2]     chrA   151-250      +
#>   [3]     chrA   301-600      +
#>   [4]     chrA     1-600      -
#>   [5]     chrA     1-600      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
cleanGaps(exampleInversion)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA   191-209      *
#>   [2]     chrA   292-300      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
