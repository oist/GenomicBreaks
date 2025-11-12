# GBreaks to xlim

Converts a `GBreaks` object to a `list` that can be passed as `xlim`
argument to
[`plotApairOfChrs`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md).

## Usage

``` r
gb2xlim(gb)
```

## Arguments

- gb:

  A `GBreaks` object.

## Value

A `list` of two pairs of numbers, where the first pair is the *start*
and *end* positions of the first
[`range`](https://rdrr.io/r/base/range.html) of the `GBreaks` object on
the *target* genome and the second pair is the *start* and *end*
positions of the same first range, on the *query* genome.

## Note

This function is intended to be run on `GBreaks` objects that have been
filtered to contain only one *seqname* on each genome.

## See also

Other genoPlotR functions:
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other Converter functions:
[`GBreaksToMatrix()`](https://oist.github.io/GenomicBreaks/reference/GBreaksToMatrix.md),
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md)

## Author

Charles Plessy

## Examples

``` r
gb2xlim(exampleColinear)
#> [[1]]
#> [1] 100 300
#> 
#> [[2]]
#> [1] 100 300
#> 
```
