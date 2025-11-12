# Index representing how the karyotype changes between chromosomes of two genomes

The index is calculated as the mean of the number of chromosomes in the
two genomes. \###############

## Usage

``` r
karyotype_index(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns a numeric value between 0 and 1. 0 is expected for same
karyotype genomes.

## Note

Here, the definition of ‘karyotype change’ is 1 - (the mean value of
change score for each chromosome in the reference genome. The change
score for chromosome Ref is calculated as 1 / (number of mapped
chromosomes in the query genome). To be noticed, the mapped length of
mapped chromosome in the query genome should be at least 10% of the
length of the reference chromosome. if chromosome Ref has only one
mapped chromosome in the query genome, the change score is 1. For 2
chromosomes, the change score is 0.5. For the symmetry, the sum of
change scores for each chromosome in query will be normalized to 1 if
the sum is more than 1.

## See also

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Examples

``` r
gb       <- GRanges(c("Ref:100-200:+",   "Ref:400-500:+",    "Re2:600-700:+"))
gb$query <- GRanges(c("Que:1100-1200:+", "Que2:1700-1800:+", "Que2:1500-1600:+"))
karyotype_index(gb)
#> [1] 0.5
```
