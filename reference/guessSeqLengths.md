# Guesstimate seqlevel lengths

When no
[`Seqinfo::seqlengths`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html)
are available, one can resort to set them as to the maximal end
coordinate found in the object.

## Usage

``` r
guessSeqLengths(gr)
```

## Arguments

- gr:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object

## Value

Returns sequence lengths that have been guessed as the maximal end
coordinates found in the `gr`, or the `gr`'s `seqlengths` if if they did
already exist.

## See also

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
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
seqlengths(exampleTranslocation$query)
#> chrB chrC 
#>  600  600 
guessSeqLengths(exampleTranslocation$query)
#> chrB chrC 
#>  600  600 
gb2 <- exampleTranslocation
seqlengths(gb2$query) <- NA
guessSeqLengths(gb2$query)
#> chrB chrC 
#>  400  300 
```
