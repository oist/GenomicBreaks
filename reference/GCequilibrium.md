# Equilibrium GC fraction (GC\*)

Computes the predicted equilibrium GC content **GC\*** from a nucleotide
substitution count matrix. GC\* is the fraction of G/C toward which base
composition will evolve under a two-class (weak vs. strong)
mutation-bias model, using the asymmetry between AT→GC (W→S) and GC→AT
(S→W) changes (Sueoka, 1962).

## Usage

``` r
GCequilibrium(m)
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

## Details

The classic directional-mutation equilibrium for GC content is:

\$\$\mathrm{GC}^\* = \frac{\text{W}\to\text{S}}{\text{W}\to\text{S}
\quad+\quad \text{S}\to\text{W}}\$\$

## Note

If the *target* genome is not the true ancestor (which is likely in
simple pairwise comparisons of extant genomes), GC\* should be
interpreted cautiously. It does not predict future GC content, but its
position relative to the GC of the *target* can still indicate the
direction of substitution bias: if GC\* is higher than the target GC,
the bias favors G/C; if lower, it favors A/T. Transpose the input matrix
to study the *query* genome the same way.

## References

Noboru Sueoka. On the genetic basis of variation and heterogeneity of
DNA base composition. (1962) *Proc Natl Acad Sci U S A* 48(4):582-92.
[doi:10.1073/pnas.48.4.582](https://doi.org/10.1073/pnas.48.4.582)

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
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
GCequilibrium(exampleSubstitutionMatrix)
#> [1] 0.5144417
GCequilibrium(t(exampleSubstitutionMatrix))
#> [1] 0.4855583
```
