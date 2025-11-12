# Convert to binned matrix of hits

Convert to binned matrix of hits

## Usage

``` r
GBreaksToMatrix(gb, ntile = 1000)
```

## Arguments

- gb:

  A `GenomicBreaks` object.

- ntile:

  The number of bins in the matrix.

## Value

A [`matrix`](https://rdrr.io/r/base/matrix.html) object.

## See also

Other Converter functions:
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md)

## Examples

``` r
m <- GBreaksToMatrix(exampleColinear5, ntile = 10)
m
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]    0    0    0    0    0    0    0    0    0     0
#>  [2,]    0    1    1    1    0    0    0    0    0     0
#>  [3,]    0    1    1    1    0    0    0    0    0     0
#>  [4,]    0    1    1    1    1    0    0    0    0     0
#>  [5,]    0    0    0    1    1    1    1    0    0     0
#>  [6,]    0    0    0    0    1    1    1    0    0     0
#>  [7,]    0    0    0    0    1    1    1    1    1     0
#>  [8,]    0    0    0    0    0    0    1    1    1     0
#>  [9,]    0    0    0    0    0    0    1    1    1     1
#> [10,]    0    0    0    0    0    0    0    0    1     1
image(m)

makeOxfordPlots(exampleColinear5)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the GenomicBreaks package.
#>   Please report the issue to the authors.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the GenomicBreaks package.
#>   Please report the issue to the authors.

```
