# Isolate the left-side gaps in inversions

For each inversion flagged by
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
returns the region of the *target* genome between the left flank of the
inversion triplet (the range carrying the flag) and its central range.

## Usage

``` r
leftInversionGaps(gb, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- detect:

  Run again
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  if `TRUE`, else reuse the `inv` flag. If the flag is absent, no gap is
  returned.

## Value

Returns a
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object representing the left-side gaps in the `GBreaks` object, one per
inversion. Its `flagIndex` metadata column contains the position in `gb`
of the left flank of each inversion.

## Note

When these two ranges are adjacent or overlap on the *target* genome,
the gap is a zero-width range that starts after the end of the left
flank.

## See also

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

## Examples

``` r
leftInversionGaps(exampleInversion)
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand | flagIndex
#>          <Rle> <IRanges>  <Rle> | <integer>
#>   [1]     chrA   191-209      * |         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Reuse flags computed with a tolerance threshold
exampleInversion |> flagInversions(tol = 19) |> leftInversionGaps(detect = FALSE)
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand | flagIndex
#>          <Rle> <IRanges>  <Rle> | <integer>
#>   [1]     chrA   191-209      * |         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleInversion |> flagInversions(tol = 18) |> leftInversionGaps(detect = FALSE)
#> GRanges object with 0 ranges and 1 metadata column:
#>    seqnames    ranges strand | flagIndex
#>       <Rle> <IRanges>  <Rle> | <integer>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
