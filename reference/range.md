# Range method for `GBreaks` objects

This is a `range` method for
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
objects, that will run
[`GenomicRanges::range`](https://rdrr.io/pkg/GenomicRanges/man/inter-range-methods.html)
on its *target* and *query* ranges and will return a new `GBreaks`
object.

## Usage

``` r
range_GBreaks(
  x,
  ...,
  with.revmap = FALSE,
  ignore.strand = FALSE,
  na.rm = FALSE
)

# S4 method for class 'GBreaks'
range(x, ..., with.revmap = FALSE, ignore.strand = FALSE, na.rm = FALSE)
```

## Arguments

- x:

  A `GBreaks` object.

- ...:

  etc

- with.revmap:

  FALSE

- ignore.strand:

  FALSE

- na.rm:

  FALSE)

## Value

tbd

## Note

`range` unconditionally ignores strand in `GBreaks` objects.

## See also

Other Bioconductor API functions:
[`GBreaks-class`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md),
[`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md),
[`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`subsetByOverlaps_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/subsetByOverlaps.md)

## Examples

``` r
range(exampleColinear3)
#> GBreaks object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-400      * | chrB:100-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
range(exampleDeletion)
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-400      * | chrB:100-300
#>   [2]     chrA   201-300      * | chrC:401-500
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
