# Set sequence lengths if there were none

Takes a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
or a
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object and if `seqlengths` are not available, sets them using the output
of the
[`guessSeqLengths`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md)
function.

## Usage

``` r
forceSeqLengths(x)

# S4 method for class 'GRanges'
forceSeqLengths(x)

# S4 method for class 'GBreaks'
forceSeqLengths(x)
```

## Arguments

- x:

  A `GBreaks` or a `GRanges` object.

## Value

`forceSeqLengths` returns the object in which the sequence lengths have
been set to the maximal end coordinates found in the object if if they
did not exist. For `GBreaks` objects it handles both the *target* and
the *query* ranges at the same time.

## See also

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

## Examples

``` r
# Prepare an example object with no seqlengths
gb <- exampleTranslocation
seqlengths(gb) <-  seqlengths(gb$query) <- NA
gb
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrC:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# Note that the new seqlengths returned by forceSeqLengths are shorter
# because we can not guess about length of the unaligned ends.
forceSeqLengths(gb)  |> seqlengths()
#> chrA 
#>  400 
exampleTranslocation |> seqlengths()
#> chrA 
#>  600 

# forceSeqLengths can take whole GBreaks objects as input, or simple GRanges
forceSeqLengths(gb)$query    |> seqlengths()
#> chrB chrC 
#>  400  300 
forceSeqLengths(gb$query)    |> seqlengths()
#> chrB chrC 
#>  400  300 
forceSeqLengths(granges(gb)) |> seqlengths()
#> chrA 
#>  400 
```
