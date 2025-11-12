# Guided scaffolding

Guided scaffolding

## Usage

``` r
scaffoldByFlipAndMerge(gr, guide, drop = FALSE)
```

## Arguments

- gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

- guide:

  A named list of data frames of contig names and orientations.

- drop:

  Drop the contigs that were not included in the guide. The new sequence
  levels will be in the same order as in the `guide`.

## See also

Other scaffolding functions:
[`flipStrandNames()`](https://oist.github.io/GenomicBreaks/reference/flipStrandNames.md),
[`longestMatchesInTarget()`](https://oist.github.io/GenomicBreaks/reference/longestMatchesInTarget.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`strandNames()`](https://oist.github.io/GenomicBreaks/reference/strandNames.md)

## Author

Michael Mansfield

Charles Plessy

## Examples

``` r
(gr <- GRanges(c("chrB:100-200:+", "chrC:201-300:+",
                 "chrB:301-400:+", "chrD:100-200:+")) |> forceSeqLengths())
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrB   100-200      +
#>   [2]     chrC   201-300      +
#>   [3]     chrB   301-400      +
#>   [4]     chrD   100-200      +
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
(g <- list(
  chrD = data.frame(contig = "chrD", orientation = 1),
  chrBC = data.frame(contig = c("chrB", "chrC"), orientation = c(1,-1))))
#> $chrD
#>   contig orientation
#> 1   chrD           1
#> 
#> $chrBC
#>   contig orientation
#> 1   chrB           1
#> 2   chrC          -1
#> 

scaffoldByFlipAndMerge(gr, g, drop = TRUE)
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chrD   100-200      +
#>   [2]    chrBC   100-200      +
#>   [3]    chrBC   301-400      +
#>   [4]    chrBC   401-500      -
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
