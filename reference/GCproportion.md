# Proportion of GC in alignment

The proportion of GC in the *target* genome in an alignment. For the
same on the *query* genome, just transpose the input matrix.

## Usage

``` r
GCproportion(m)
```

## Arguments

- m:

  A matrix of counts or probabilities for bases of the *target* genome
  to be aligned to bases on the *query* genome. As a convenience it can
  also receive a list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

## Value

A numeric value between 0 and 1.

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
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
GCproportion(exampleSubstitutionMatrix)
#> [1] 0.4347341
GCproportion(t(exampleSubstitutionMatrix))
#> [1] 0.4302136
```
