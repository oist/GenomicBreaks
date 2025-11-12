# Find longest matches from target to query genome.

Using a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object representing the alignment of a *query* genome on a *target*
genome, finds the longest match of each sequence level (representing
contigs, scaffolds, etc.) of the *query* on the *target*.

## Usage

``` r
longestMatchesInTarget(gb, min.width = 10000, min.matches = 2)
```

## Arguments

- gb:

  A `GBreaks` object

- min.width:

  Minimum width of a match (on the query genome) for being considered.

- min.matches:

  Discard query sequences that have fewer longest matches than
  `min.matches` on the target. Default is 2, so that only results
  relevant to chaining genomes are kept.

## Value

Returns a
[`GenomicRanges::GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html)
object containing one
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object per sequence on the query genome.

## Details

Each sequence of the *query* is represented only once in the output, but
sequences of the *target* genome can be represented multiple times if
they are the longest match of multiple *query* genome sequences. When
the *target* genome is more contiguous than the *query* genome, and if
there are no major structural variations between them, this will reveal
arrangements of colinear sequences in the query genome.

For a more compact version of the results, the output of this function
can be piped to `strandNames(query = TRUE)`.

## See also

Other scaffolding functions:
[`flipStrandNames()`](https://oist.github.io/GenomicBreaks/reference/flipStrandNames.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`scaffoldByFlipAndMerge()`](https://oist.github.io/GenomicBreaks/reference/scaffoldByFlipAndMerge.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`strandNames()`](https://oist.github.io/GenomicBreaks/reference/strandNames.md)

## Author

Charles Plessy

## Examples

``` r
exampleColinear3
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrB:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleColinear3 |> longestMatchesInTarget(min.width = 0, min.matches = 1)
#> GRangesList object of length 1:
#> $chrB
#> GBreaks object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrB   100-200      + | chrA:100-200
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
#> 
```
