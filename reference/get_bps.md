# Get breakpoints

Given a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
or
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object, the function produces a `GRanges` object detailing the
breakpoints only.

## Usage

``` r
get_bps(
  gr,
  direction = c("both", "left", "right", "mid"),
  stranded = FALSE,
  sorted = TRUE
)
```

## Arguments

- gr:

  `GRanges` object containing pairwise alignment

- direction:

  Return the breakpoints on `both`, `left` or `right` side(s) of the
  range, or at the `mid`point between ranges.

- stranded:

  If `TRUE`, will assign a `+` strand to the left-side breakpoints and a
  `-` strand to the right-side ones.

- sorted:

  Sorts the object before returning it.

## Value

`GRanges` object of the breakpoints

## See also

Other Reducing functions:
[`cleanGaps()`](https://oist.github.io/GenomicBreaks/reference/cleanGaps.md),
[`wholeGenomeClassification()`](https://oist.github.io/GenomicBreaks/reference/wholeGenomeClassification.md)

## Examples

``` r
get_bps(exampleInversion)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA       100      *
#>   [2]     chrA       190      *
#>   [3]     chrA       210      *
#>   [4]     chrA       291      *
#>   [5]     chrA       301      *
#>   [6]     chrA       400      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
get_bps(exampleInversion, direction = "left")
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA       100      *
#>   [2]     chrA       210      *
#>   [3]     chrA       301      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
get_bps(exampleInversion, stranded = TRUE)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA       100      +
#>   [2]     chrA       190      -
#>   [3]     chrA       210      +
#>   [4]     chrA       291      -
#>   [5]     chrA       301      +
#>   [6]     chrA       400      -
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
get_bps(exampleInversion, direction = "right", stranded = TRUE)
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA       190      -
#>   [2]     chrA       291      -
#>   [3]     chrA       400      -
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
get_bps(exampleInversion, direction = "mid")
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA       200      *
#>   [2]     chrA       296      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
