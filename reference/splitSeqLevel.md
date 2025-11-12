# Split a seqlevel in two pieces

For scaffolding or plotting purposes, it may be useful to split some
sequences into smaller ones.

## Usage

``` r
splitSeqLevel(gb, seq, bp)
```

## Arguments

- gb:

  A `GBreaks` object.

- seq:

  A one seqlevel from `gr`

- bp:

  The breakpoint where to split.

## Value

Returns a modified `GRanges` object in which the sequence has been split
Its [`Seqinfo::seqinfo`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html)
has a new entry for the new levels, and the old level is not removed. If
no `seqlengths` were present in the original object, they are
arbitrarily set as the maximal end value for each `seqlevel`.

## Note

This function only splits at breakpoints.

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
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

Other scaffolding functions:
[`flipStrandNames()`](https://oist.github.io/GenomicBreaks/reference/flipStrandNames.md),
[`longestMatchesInTarget()`](https://oist.github.io/GenomicBreaks/reference/longestMatchesInTarget.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`scaffoldByFlipAndMerge()`](https://oist.github.io/GenomicBreaks/reference/scaffoldByFlipAndMerge.md),
[`strandNames()`](https://oist.github.io/GenomicBreaks/reference/strandNames.md)

## Author

Charles Plessy

## Examples

``` r
splitSeqLevel(exampleInsertion, "chrA", 200)
#> GBreaks object with 3 ranges and 1 metadata column:
#>               seqnames    ranges strand |        query
#>                  <Rle> <IRanges>  <Rle> |    <GRanges>
#>          chrC     chrC   401-500      + | chrB:201-300
#>   chrA.chrA_1   chrA_1   100-200      + | chrB:100-200
#>   chrA.chrA_2   chrA_2     1-100      + | chrB:301-400
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome
```
