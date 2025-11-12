# Proportion of gaps in alignment

The proportion of nucleotides that are part of an alignment but not
matched to another nucleotide. If you divide it by the percent
difference, you obtain the S/I ratio studied by Chen and coll., 2009.

## Usage

``` r
gapProportion(m)
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

## Note

This is not a measure of distance. At greater evolutionary divergence,
genome alignment becomes increasingly difficult, leaving only the most
conserved regions alignable. Consequently, the proportion of gaps tends
to decrease.

## References

Chen JQ, Wu Y, Yang H, Bergelson J, Kreitman M, Tian D. (2009).
Variation in the ratio of nucleotide substitution and indel rates across
genomes in mammals and bacteria. *Mol Biol Evol.* 26(7):1523-31.
[doi:10.1093/molbev/msp063](https://doi.org/10.1093/molbev/msp063)

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md)

## Author

Charles Plessy

## Examples

``` r
gapProportion(exampleSubstitutionMatrix)
#> [1] 0.1397383
```
