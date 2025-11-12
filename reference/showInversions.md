# Show inversions and their flanking blocks.

Show inversions and their flanking blocks.

## Usage

``` r
showInversions(gb, rename = TRUE, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md).

- rename:

  Replace range names by their numeric order before subsetting.

- detect:

  Run again
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  if `TRUE`, else reuse the `inv` flag or fail if it is absent.

## Value

Returns the `GBreaks` object in which all ranges that are not part of an
inversion triplet have been discarded. If the object was missing the
`inv` metadata column, return the object after discarding all of its
ranges.

## See also

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md)

## Examples

``` r
showInversions(exampleInversion)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |        query   inv
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   1     chrA   100-190      + | chrB:100-190  TRUE
#>   2     chrA   210-291      - | chrB:210-291 FALSE
#>   3     chrA   301-400      + | chrB:301-400 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
