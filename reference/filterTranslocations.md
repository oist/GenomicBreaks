# Extract central blocks in translocations

Extract central blocks in translocations

## Usage

``` r
filterTranslocations(gb, rename = TRUE, remove = FALSE, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- rename:

  Replace range names by their numeric order before subsetting.

- remove:

  Filter out instead of filtering in.

- detect:

  Run again
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  if `TRUE`, else reuse the `tra` flag or fail if it is absent.

## Value

Returns the `GBreaks` object in which all ranges that are not the
central part of an inversion triplet have been discarded. If the object
was missing the `tra` metadata column, return the object after
discarding all of its ranges.

## See also

Other Translocation functions:
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md)

## Examples

``` r
filterTranslocations(exampleTranslocation)
#> GBreaks object with 1 range and 2 metadata columns:
#>     seqnames    ranges strand |        query   tra
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   2     chrA   201-300      + | chrC:201-300 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
filterTranslocations(remove = TRUE, exampleTranslocation)
#> GBreaks object with 2 ranges and 2 metadata columns:
#>     seqnames    ranges strand |        query   tra
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   1     chrA   100-200      + | chrB:100-200  TRUE
#>   3     chrA   301-400      + | chrB:301-400 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
