# Reverse genomic ranges

Reverse genomic ranges by flipping the strand and moving the origin of
the coordinate system to the opposite side.

## Usage

``` r
reverse(x, ...)

# S4 method for class 'GRanges'
reverse(x, ...)

# S4 method for class 'GBreaks'
reverse(x, query = FALSE, ...)
```

## Arguments

- x:

  A `GBreaks` or a `GRanges` object

- ...:

  Additional arguments (ignored)

- query:

  On `GBreaks` objects, operate on on the *query* or the *target*.

## Value

Returns the modified `GBreaks` or `GRanges` object.

## See also

See also the
[`IRanges::reverse`](https://rdrr.io/pkg/IRanges/man/reverse-methods.html)
function.

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

Other Bioconductor API functions:
[`GBreaks-class`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md),
[`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md),
[`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md),
[`range_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/range.md),
[`subsetByOverlaps_GBreaks()`](https://oist.github.io/GenomicBreaks/reference/subsetByOverlaps.md)

## Author

Charles Plessy

## Examples

``` r
exampleInsertion
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrB:301-400
#>   [3]     chrC   401-500      + | chrB:201-300
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
reverse(exampleInsertion)
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   401-501      - | chrB:100-200
#>   [2]     chrA   301-400      - | chrB:301-400
#>   [3]     chrC   101-200      - | chrB:201-300
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
reverse(exampleInsertion, query = TRUE)
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      - | chrB:401-501
#>   [2]     chrA   201-300      - | chrB:201-300
#>   [3]     chrC   401-500      - | chrB:301-400
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
