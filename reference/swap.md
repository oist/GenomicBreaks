# Swap target and query genomes

Produce a new object in which the information about the *target* and
*query* genomes have been inverted. As a large number of *GenomicBreaks*
functions require `GBreaks` objects that were sorted with strand
information ignored, a `sort` option is provided to do that.

## Usage

``` r
swap(gb, sort = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- sort:

  Sort the returned object, ignoring strand.

## Value

A new `GBreaks` object representing the input's *query* ranges, with the
input's *target* ranges in the *query* slot. The strand information is
unchanged.

## See also

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
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md)

## Examples

``` r
swap(exampleDeletion)
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrB   100-200      + | chrA:100-200
#>   [2]     chrC   401-500      + | chrA:201-300
#>   [3]     chrB   201-300      + | chrA:301-400
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
swap(exampleDeletion, sort = TRUE)
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrB   100-200      + | chrA:100-200
#>   [2]     chrB   201-300      + | chrA:301-400
#>   [3]     chrC   401-500      + | chrA:201-300
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
