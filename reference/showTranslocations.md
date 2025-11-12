# Show translocations and their flanking blocks.

Show translocations and their flanking blocks.

## Usage

``` r
showTranslocations(gb, rename = TRUE, detect = TRUE)
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
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  if `TRUE`, else reuse the `tra` flag or fail if it is absent.

## Value

Returns the `GBreaks` object in which all ranges that are not part of a
translocation triplet have been discarded. If the object was missing the
`tra` metadata column, return the object after discarding all of its
ranges.

## See also

Other Translocation functions:
[`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md)

## Examples

``` r
showTranslocations(exampleTranslocation)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |        query   tra
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   1     chrA   100-200      + | chrB:100-200  TRUE
#>   2     chrA   201-300      + | chrC:201-300 FALSE
#>   3     chrA   301-400      + | chrB:301-400 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
