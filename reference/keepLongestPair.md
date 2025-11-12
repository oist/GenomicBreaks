# Keep the longest pair of features

This feature is useful when analysing genomes comprising a main
chromosomes and extra plasmids or smaller chromosomes. Assuming that the
main chroomosome can be detected only by its length, this function
extracts the alignment between the pair of main chromosomes only.

## Usage

``` r
keepLongestPair(gb, drop = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object

- drop:

  Drop unused sequence levels.

## Value

A `GBreaks` object with only one sequence feature in use per genome.

## See also

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

## Author

Charles Plessy

## Examples

``` r
exampleTranslocation
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrC:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleTranslocation |> keepLongestPair()
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
