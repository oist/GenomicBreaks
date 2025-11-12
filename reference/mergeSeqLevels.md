# Merge seqlevels in a larger one

For scaffolding or plotting purposes, it may be useful to merge some
sequences into larger ones.

## Usage

``` r
mergeSeqLevels(gr, seqs, name)

mergeSeqLevels_to_DF(gr, seqs, name)
```

## Arguments

- gr:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

- seqs:

  A character vector of
  [`Seqinfo::seqlevels`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html)
  from `gr`

- name:

  The name of the new sequence level to be added

## Value

Returns a modified `GRanges` object in which the sequences have been
merged. Its
[`Seqinfo::seqinfo`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html) has a
new entry for the new level, and the old levels are not removed. If no
`seqlengths` were present in the original object, they are arbitrarily
set as the maximal end value for each `seqlevel`.

The `mergeSeqLevels_to_DF` function returns a `DataFrame` in which the
`start` and `end` columns are in `numeric` mode. This is to cirvumvent
the fact that `GenomicRanges` object hardcode the mode of *start* and
*end* positions to `integer`, which does not allow values larger than
2,147,483,647, which does not allow to merge sequence levels of
mammalian or larger-scale genomes.

## Note

Be careful that in some cases it is needed to "flip" the sequence
feature with
[`reverse`](https://oist.github.io/GenomicBreaks/reference/reverse.md)
before merging, for instance when colinearity is with its reverse
strand.

## See also

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

Other scaffolding functions:
[`flipStrandNames()`](https://oist.github.io/GenomicBreaks/reference/flipStrandNames.md),
[`longestMatchesInTarget()`](https://oist.github.io/GenomicBreaks/reference/longestMatchesInTarget.md),
[`scaffoldByFlipAndMerge()`](https://oist.github.io/GenomicBreaks/reference/scaffoldByFlipAndMerge.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`strandNames()`](https://oist.github.io/GenomicBreaks/reference/strandNames.md)

## Author

Charles Plessy

## Examples

``` r
gb       <- GRanges(c("XSR:101-180:+", "XSR:201-300:+",  "XSR:320-400:+"))
gb$query <- GRanges(c( "S1:101-200",      "S2:1-100",    "S3:1-100"))
seqlengths(gb$query) <- c(200, 100, 100)
genome(gb$query) <- "GenomeX"
isCircular(gb$query) <- rep(FALSE, 3)
seqinfo(gb$query)
#> Seqinfo object with 3 sequences from GenomeX genome:
#>   seqnames seqlengths isCircular  genome
#>   S1              200      FALSE GenomeX
#>   S2              100      FALSE GenomeX
#>   S3              100      FALSE GenomeX
gb <- GBreaks(gb)
gb$query <- mergeSeqLevels(gb$query, c("S2", "S3"), "Scaf1")
gb
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |         query
#>          <Rle> <IRanges>  <Rle> |     <GRanges>
#>   [1]      XSR   101-180      + |    S1:101-200
#>   [2]      XSR   201-300      + |   Scaf1:1-100
#>   [3]      XSR   320-400      + | Scaf1:101-200
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
seqinfo(gb$query)
#> Seqinfo object with 4 sequences from GenomeX genome:
#>   seqnames seqlengths isCircular  genome
#>   Scaf1           200      FALSE GenomeX
#>   S1              200      FALSE GenomeX
#>   S2              100      FALSE GenomeX
#>   S3              100      FALSE GenomeX

mergeSeqLevels(gb, seqlevelsInUse(gb), "AllMerged")
#> GBreaks object with 3 ranges and 1 metadata column:
#>        seqnames    ranges strand |         query
#>           <Rle> <IRanges>  <Rle> |     <GRanges>
#>   [1] AllMerged   101-180      + |    S1:101-200
#>   [2] AllMerged   201-300      + |   Scaf1:1-100
#>   [3] AllMerged   320-400      + | Scaf1:101-200
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
