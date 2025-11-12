# Convert `GRanges` object to `dna_seg` format

The
[`genoPlotR::dna_seg`](https://rdrr.io/pkg/genoPlotR/man/dna_seg.html)
class represents *DNA segments* using a `data.frame` format.

## Usage

``` r
gr2dna_seg(gr, ...)
```

## Arguments

- gr:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

- ...:

  Extra arguments passed to `dna_seg()`.

## Value

Returns a `genoPlotR::dna_seq` or `NULL` if the input object has a
length of zero.

## See also

Other genoPlotR functions:
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other Converter functions:
[`GBreaksToMatrix()`](https://oist.github.io/GenomicBreaks/reference/GBreaksToMatrix.md),
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md)

## Author

Charles Plessy

## Examples

``` r
gr2dna_seg(exampleInversion)
#>             name start end strand  col fill lty lwd pch cex gene_type
#> 1 chrA:100-190:+   100 190      1 blue blue   1   1   8   1    blocks
#> 2 chrA:210-291:-   210 291     -1 blue blue   1   1   8   1    blocks
#> 3 chrA:301-400:+   301 400      1 blue blue   1   1   8   1    blocks
gr2dna_seg(exampleInversion$query)
#>           name start end strand  col fill lty lwd pch cex gene_type
#> 1 chrB:100-190   100 190      1 blue blue   1   1   8   1    blocks
#> 2 chrB:210-291   210 291      1 blue blue   1   1   8   1    blocks
#> 3 chrB:301-400   301 400      1 blue blue   1   1   8   1    blocks
```
