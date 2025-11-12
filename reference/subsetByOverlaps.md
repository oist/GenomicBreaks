# subsetByOverlaps method for `GBreaks` objects

This is a `subsetByOverlaps` method for
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
objects, that will run
[`GenomicRanges::subsetByOverlaps`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)
on its *target* and *query* ranges and will return a new `GBreaks`
object.

## Usage

``` r
subsetByOverlaps_GBreaks(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...
)

# S4 method for class 'GBreaks,GBreaks'
subsetByOverlaps(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...
)
```

## Arguments

- x:

  A `GBreaks` object.

- ranges:

  Another `GBreaks` object.

- maxgap:

  tbd

- minoverlap:

  tbd

- type:

  tbd

- invert:

  tbd

- ...:

  etc

## Value

tbd

## Note

`range` unconditionally ignores strand in `GBreaks` objects.

## See also

Other Bioconductor API functions:
[`GBreaks-class`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md),
[`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md),
[`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md),
[`range_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/range.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md)

## Examples

``` r
subsetByOverlaps(exampleColinear3, exampleColinear3)
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrB:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
