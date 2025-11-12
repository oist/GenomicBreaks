# GC pressure

Computes the following equation: \$\$\frac{\text{W}\to\text{S}
\quad-\quad \text{S}\to\text{W}}{\text{W}\to\text{S} \quad+\quad
\text{S}\to\text{W}}\$\$

## Usage

``` r
GCpressure(m)
```

## Arguments

- m:

  A matrix of counts or probabilities for bases of the *target* genome
  to be aligned to bases on the *query* genome. As a convenience it can
  also receive a list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

## Value

A numeric value between \\-1\\ and \\1\\.

## Note

This index was suggested to by by M365 Copilot while using it to write
the
[`GCequilibrium`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md)
function and its documentation. Despite my attempts I could not find a
paper that uses it explicitely… I am including it in this package
because I am screening as many indices that I can extract from an
alignment pair.

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md)

## Author

Charles Plessy

## Examples

``` r
GCpressure(exampleSubstitutionMatrix)
#> [1] 0.02888339
GCpressure(t(exampleSubstitutionMatrix))
#> [1] -0.02888339
```
