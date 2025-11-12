# Isolate the left-side gaps in inversions

Isolate the left-side gaps in inversions

## Usage

``` r
leftInversionGaps(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns a
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object representing the left-side gaps in the `GBreaks` object.

## See also

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

## Examples

``` r
leftInversionGaps(flagInversions(exampleInversion))
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrA   191-209      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
