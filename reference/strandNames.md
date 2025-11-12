# Sequence name with strand information

Extract sequence names and pastes strand information to it.

## Usage

``` r
strandNames(gb, flip = FALSE, query = FALSE)
```

## Arguments

- gb:

  A `GBreaks` object or a
  [`GenomicRanges::GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html)
  of `GBreaks` objects.

- flip:

  Flip the strand names.

- query:

  Work on the query genome.

## Value

Returns a character vector, or a list of character vectors if the input
was a
[`GenomicRanges::GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html).

## See also

Other scaffolding functions:
[`flipStrandNames()`](https://oist.github.io/GenomicBreaks/reference/flipStrandNames.md),
[`longestMatchesInTarget()`](https://oist.github.io/GenomicBreaks/reference/longestMatchesInTarget.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`scaffoldByFlipAndMerge()`](https://oist.github.io/GenomicBreaks/reference/scaffoldByFlipAndMerge.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md)

## Author

Charles Plessy

## Examples

``` r
strandNames(exampleColinear)
#> [1] "chrA+" "chrA+"
strandNames(exampleColinear, query = TRUE)
#> [1] "chrB+" "chrB+"
```
