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
  if `TRUE`, else reuse the `tra` and `trs` flags. If they are absent,
  an empty object is returned.

## Value

Returns the `GBreaks` object in which all ranges that are not part of a
payload region (`U A V`), of a site (`X Y`), or the moved range of a
site, have been discarded (see
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)).

## See also

Other Translocation functions:
[`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md)

## Examples

``` r
showTranslocations(exampleTranslocation)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>     seqnames    ranges strand |        query   tra   trs
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   1     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   2     chrA   201-300      + | chrC:201-300 FALSE FALSE
#>   3     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
showTranslocations(swap(exampleTranslocation, sort = TRUE))
#> GBreaks object with 3 ranges and 3 metadata columns:
#>     seqnames    ranges strand |        query   tra   trs
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   1     chrB   100-200      + | chrA:100-200 FALSE  TRUE
#>   2     chrB   301-400      + | chrA:301-400 FALSE FALSE
#>   3     chrC   201-300      + | chrA:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
