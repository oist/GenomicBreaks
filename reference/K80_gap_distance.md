# Kimura's 2-parameter distance extended with gaps

Nishimaki and Sato's K2P + Gap distance extends Kimura's 2-parameter
distance by treating gaps as insertion/deletion states instead of
removing gap-containing sites.

## Usage

``` r
K80_gap_distance(m)
```

## Arguments

- m:

  A matrix of counts or probabilities for bases of the *target* genome
  to be aligned to bases on the *query* genome. As a convenience it can
  also receive a list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

## Value

Returns a numeric value showing the evolutionary distance between two
genomes. The larger the value, the more different the two genomes are.

## Details

The distance is calculated as \$\$ K = \frac{3}{4} w \log(w) -
\frac{w}{2} \log\left((S - P) \sqrt{S + P - Q}\right) \$\$ where \\S\\
is the probability of identical nucleotide pairs, \\P\\ is the
probability of transition-type nucleotide pairs, \\Q\\ is the
probability of transversion-type nucleotide pairs, and \\w\\ is the
nucleotide occupancy probability.

When there are no gaps, \\w = 1\\, and this expression reduces to the
usual K80 distance.

## References

Kimura, M. (1980). "A simple method for estimating evolutionary rates of
base substitutions through comparative studies of nucleotide sequences."
*Journal of Molecular Evolution*, 16, 111–120. DOI:
[doi:10.1007/BF01731581](https://doi.org/10.1007/BF01731581)

Nishimaki, T. and Sato, K. (2019). "An Extension of the Kimura
Two-Parameter Model to the Natural Evolutionary Process." *Journal of
Molecular Evolution*, 87, 60–67. DOI:
[doi:10.1007/s00239-018-9885-1](https://doi.org/10.1007/s00239-018-9885-1)

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`JC69_distance_allseq()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance_allseq.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`JC69_distance_allseq()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance_allseq.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Charles Plessy

## Examples

``` r
K80_gap_distance(exampleSubstitutionMatrix)
#> [1] 0.313952

# When there are no gaps, it returns the same as the K80 distance
nogaps <- exampleSubstitutionMatrix
nogaps["-",] <- 0
nogaps[,"-"] <- 0
K80_gap_distance(nogaps)
#> [1] 0.2789688
K80_distance(nogaps)
#> [1] 0.2789688
```
