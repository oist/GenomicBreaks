# Extract moved blocks in translocations

Extract moved blocks in translocations

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
  if `TRUE`, else reuse the `tra` and `trs` flags. If they are absent,
  an empty object is returned.

## Value

Returns the `GBreaks` object in which all ranges that are not moved
ranges have been discarded: the centre of each payload region, and the
range found between the two ranges of each site on the *query* genome
(see
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)).
A range found by both flags is returned once. With `remove = TRUE`, the
moved ranges are discarded instead.

## See also

Other Translocation functions:
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md)

## Examples

``` r
filterTranslocations(exampleTranslocation)
#> GBreaks object with 1 range and 3 metadata columns:
#>     seqnames    ranges strand |        query   tra   trs
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   2     chrA   201-300      + | chrC:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
filterTranslocations(remove = TRUE, exampleTranslocation)
#> GBreaks object with 2 ranges and 3 metadata columns:
#>     seqnames    ranges strand |        query   tra   trs
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   1     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   3     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
filterTranslocations(swap(exampleTranslocation, sort = TRUE))
#> GBreaks object with 1 range and 3 metadata columns:
#>     seqnames    ranges strand |        query   tra   trs
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   3     chrC   201-300      + | chrA:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
