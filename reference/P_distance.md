# Percent difference Distance

The P-distance is simply the percentage of the aligned base pairs where
the two bases are different. See Details for a discussion on the many
ways to compute percent identity (and therefore difference.)

The `P_error` function reports the standard error given the alignment
length \\n\\ and percent difference \\p\\ computed from the matrix:
\\\sqrt{p(1 - p) / n}\\.

## Usage

``` r
P_distance(m, denominator = c("L3", "L1", "L2", "L4"))

P_error(m, denominator = c("L3", "L1", "L2", "L4"))
```

## Arguments

- m:

  A matrix of counts (`P_distance`, `P_variance`) or probabilities
  (`P_distance` only) for bases of the *target* genome to be aligned to
  bases on the *query* genome. As a convenience it can also receive a
  list produced by the
  [`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
  function, containing this matrix.

- denominator:

  Denominator according to the nomenclature of May (2004). Default is
  `L3`.

## Value

A numeric value representing the evolutionary distance between two
genomes. The larger the value, the more divergent the genomes.

## Details

The P-distance is computed as 1 minus the percent sequence identity
according to the definitions of May (2004), which is the ratio between
the number of matching bases and four possible denominators in use at
that time.

- L1: Length of the shorter sequence.

- L2: Number of aligned positions, i.e., alignment length (includes
  gaps, if any).

- L3: Number of aligned residue pairs, i.e., identities and
  nonidentities (excludes gaps, if any).

- L4: Arithmetic mean sequence length.

## References

May, A. C. W. (2004). Percent Sequence Identity: The Need to Be
Explicit. *Structure*, 12(5), 737–738.
[doi:10.1016/j.str.2004.04.001](https://doi.org/10.1016/j.str.2004.04.001)

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
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
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md)

## Author

Charles Plessy

## Examples

``` r
P_distance(exampleSubstitutionMatrix)
#> [1] 0.2312837

# example code
P_error(exampleSubstitutionMatrix)
#> [1] 7.701531e-05
```
