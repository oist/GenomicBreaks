# Gene Order Conservation

The Gene Order Conservation (GOC) number defined by Rocha (2003) is:
“*the average number of orthologues for which the consecutive orthologue
co-occurs close by in the other genome. It varies between 0 (no
co-occurrence) and 1 (complete gene order conservation)*”.

## Usage

``` r
GOC(gb, vicinity = 4, debug = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- vicinity:

  How far to search for the neighbour orthologue.

- debug:

  See below.

## Value

Returns a numeric value between 0 and 1. If `debug = TRUE`, returns a
copy of the `gb` object with additional columns showing details of the
computation.

## Note

Note that calculating GOC on whole-genome alignments is not expected to
produce meaningful results. This function is more useful when comparing
the position of orthologues, represented in a `GBreaks` object.

## References

Rocha, Eduardo P C. “DNA repeats lead to the accelerated loss of gene
order in bacteria.” *Trends in genetics : TIG* vol. 19,11 (2003): 600-3.
[doi:10.1016/j.tig.2003.09.011](https://doi.org/10.1016/j.tig.2003.09.011)

## See also

Other Colinearity functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Examples

``` r
exampleColinear
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-150      + | chrB:100-150
#>   [2]     chrA   251-300      + | chrB:251-300
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
GOC(exampleColinear)
#> [1] 1

exampleTranslocation
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrC:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
GOC(exampleTranslocation)
#> [1] 0.5
GOC(exampleTranslocation, v=1)
#> [1] 0

# GOC computation is strandless
exampleInversion
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   210-291      - | chrB:210-291
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
GOC(exampleInversion)
#> [1] 1
```
