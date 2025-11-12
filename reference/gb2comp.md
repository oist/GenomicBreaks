# Convert `GBreaks` object to a `comparison` object

Convert `GBreaks` object to a `comparison` object

## Usage

``` r
gb2comp(gb, color = NULL, ignore.strand = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- color:

  (optional) A vector of colors, of same length as `gb`.

- ignore.strand:

  Ignore strand information?

## Value

Returns a
[`genoPlotR::comparison`](https://rdrr.io/pkg/genoPlotR/man/comparison.html)
object.

## Note

Note that the `comparison` objects assume that the ranges are on the
same sequence feature. It is your responsibility that the input
`GBreaks` object also conforms to this assumption.

## See also

Other genoPlotR functions:
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other Converter functions:
[`GBreaksToMatrix()`](https://oist.github.io/GenomicBreaks/reference/GBreaksToMatrix.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md)

## Author

Charles Plessy

## Examples

``` r
gb2comp(exampleInversion)
#>   start1 end1 start2 end2 direction      color            col
#> 1    100  190    100  190         1 darksalmon     darksalmon
#> 2    210  291    291  210        -1 darksalmon cornflowerblue
#> 3    301  400    301  400         1 darksalmon     darksalmon
```
