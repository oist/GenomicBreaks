# Tamura 1992 distance

The Tamura 1992 (T92) distance extends the K80 distance by taking `GC`
content into account. It is calculated as \\-h \ln \left(1 -
\frac{p}{h} - q\right) - \frac{1}{2} \times (1 - h) \ln\left(1 - 2
q\right)\\, where \\p\\ is the probability of transition, \\q\\ the
probability of transversion, \\h = 2\theta (1 - \theta)\\ and \\\theta\\
is the `GC` content. See the
[Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#T92_model_(Tamura_1992))
for more details.

## Usage

``` r
T92_distance(m, gc = c("average", "target", "query"))
```

## Arguments

- m:

  A matrix of counts or probabilities for bases of the *target* genome
  to be aligned to bases on the *query* genome. As a convenience it can
  also receive a list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

- gc:

  Calculate the GC content from the *target*, the *query* or average
  both?

## Value

Returns a numeric value show the evolutionary distance between two
genomes. the larger the value, the more different the two genomes are.

## References

Tamura, K. (1992). "Estimation of the number of nucleotide substitutions
when there are strong transition-transversion and G+C-content biases."
*Molecular Biology and Evolution*, 9(4), 678–687. DOI:
[10.1093/oxfordjournals.molbev.a040752](https://doi.org/10.1093/oxfordjournals.molbev.a040752)

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
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Zikun Yang

Charles Plessy

## Examples

``` r
T92_distance(exampleSubstitutionMatrix)
#> [1] 0.2794627
T92_distance(exampleSubstitutionMatrix, gc="target")
#> [1] 0.2795033
T92_distance(exampleSubstitutionMatrix, gc="query")
#> [1] 0.2794238
```
