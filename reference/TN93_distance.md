# Tamura–Nei 1993 distance

The Tamura–Nei 1993 (TN93) distance extends the T92 distance by allowing
unequal base frequencies and two distinct transition components, plus
transversions.

## Usage

``` r
TN93_distance(m)
```

## Arguments

- m:

  A matrix of counts or probabilities for bases of the *target* genome
  aligned to bases on the *query* genome. As a convenience it can also
  receive a list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

## Value

Returns a numeric value show the evolutionary distance between two
genomes. the larger the value, the more different the two genomes are.

## Details

Let \\P_R\\ be the fraction of purine transitions (A↔G), \\P_Y\\ the
fraction of pyrimidine transitions (C↔T), and \\Q\\ the fraction of
transversions among A/C/G/T pairs. Let \\\pi_A, \pi_C, \pi_G, \pi_T\\ be
the average base frequencies across the two sequences (row and column
marginals averaged), \\\pi_R = \pi_A + \pi_G\\, and \\\pi_Y = \pi_C +
\pi_T\\.

The TN93 distance is: \$\$ d\_{\mathrm{TN93}} = -\\2 \pi_A \pi_G \\
\ln\\\left(1 - \frac{P_R}{2 \pi_A \pi_G}\right) -\\2 \pi_C \pi_T \\
\ln\\\left(1 - \frac{P_Y}{2 \pi_C \pi_T}\right) -\\2 \pi_R \pi_Y \\
\ln\\\left(1 - \frac{Q}{2 \pi_R \pi_Y}\right) \$\$ v

## References

Tamura, K. and Nei, M. (1993). Estimation of the number of nucleotide
substitutions in the control region of mitochondrial DNA in humans and
chimpanzees. *Molecular Biology and Evolution* **10**:512–26.
[doi:10.1093/oxfordjournals.molbev.a040023](https://doi.org/10.1093/oxfordjournals.molbev.a040023)

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
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Charles Plessy / M365 Copilot (GPT-5 on)

## Examples

``` r
TN93_distance(exampleSubstitutionMatrix)
#> [1] 0.288648
```
