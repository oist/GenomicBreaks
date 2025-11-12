# Extract central blocks in double inversions

Extract central blocks in double inversions

## Usage

``` r
filterDoubleInversions(gb, rename = TRUE, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- rename:

  Replace range names by their numeric order before subsetting.

- detect:

  Run again
  [`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md)
  if `TRUE`, else reuse the `Dbl` flag or fail if it is absent.

## Value

Returns the `GBreaks` object in which all ranges that are not the
central part of a double inversion quintuplet have been discarded. If
the object was missing the `Dbl` metadata column, return the object
after discarding all of its ranges.

## See also

Other Inversion functions:
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

## Examples

``` r
filterDoubleInversions(exampleDoubleInversion1)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |        query       Dbl
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   3     chrA   300-390      - | chrB:200-290     FALSE
#>   4     chrA   400-490      - | chrB:300-390     FALSE
#>   5     chrA   500-590      + | chrB:500-590     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
