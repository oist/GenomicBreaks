# Convert `GBreaks` object to a list of `dna_seg` objects

Convert `GBreaks` object to a list of `dna_seg` objects

## Usage

``` r
gb2dna_seg(gb, ...)
```

## Arguments

- gb:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  or a
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- ...:

  Extra arguments passed to `dna_seg()`.

## Value

Returns a `list` of two
[`genoPlotR::dna_seg`](https://rdrr.io/pkg/genoPlotR/man/dna_seg.html)
objects, respectively for the *target* ranges the *query* ranges.

## See also

Other genoPlotR functions:
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other Converter functions:
[`GBreaksToMatrix()`](https://oist.github.io/GenomicBreaks/reference/GBreaksToMatrix.md),
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md)

## Author

Charles Plessy

## Examples

``` r
gb2dna_seg(exampleInversion)
#> $target
#>             name start end strand  col fill lty lwd pch cex gene_type
#> 1 chrA:100-190:+   100 190      1 blue blue   1   1   8   1    blocks
#> 2 chrA:210-291:-   210 291      1 blue blue   1   1   8   1    blocks
#> 3 chrA:301-400:+   301 400      1 blue blue   1   1   8   1    blocks
#> 
#> $query
#>           name start end strand  col fill lty lwd pch cex gene_type
#> 1 chrB:100-190   100 190      1 blue blue   1   1   8   1    blocks
#> 2 chrB:210-291   210 291     -1 blue blue   1   1   8   1    blocks
#> 3 chrB:301-400   301 400      1 blue blue   1   1   8   1    blocks
#> 
```
