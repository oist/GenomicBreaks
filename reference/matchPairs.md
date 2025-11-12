# Match pairs of sequences

For each sequence of the *target* genome, assign a sequence of the
*query* genome, and discard all the ranges that are not from matched
pairs.

## Usage

``` r
matchPairs(gb, drop = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- drop:

  Drop unused sequence levels.

## Value

A `GBreaks` object.

## Note

This is not a reciprocal best match function, because it assumes that
the genomes are close enough and that there were no major karyotype
changes, so that the pairing is a trivial operation.

## See also

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
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
matchPairs(exampleTranslocation)
#> GBreaks object with 2 ranges and 1 metadata column:
#>        seqnames    ranges strand |        query
#>           <Rle> <IRanges>  <Rle> |    <GRanges>
#>   chrA     chrA   100-200      + | chrB:100-200
#>   chrA     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
