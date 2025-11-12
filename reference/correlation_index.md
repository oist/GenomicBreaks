# Index representing how correlated syntenic regions are

A sequence feature on the *target* genome is considered ‘syntenic’ to
the feature of the *query* genome matched by the most of its alignments.
The correlation coefficient of the center of syntenic alignments is then
calculated. Then, the computed numbers for all the *target*'s features
are averaged as a mean weighted by is the total number of aligned bases
of each feature. Thus, a number close to 1 is expected for closely
related genomes.

## Usage

``` r
correlation_index(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns a numeric value between 0 and 1.

## Note

Here, the definition of ‘synteny’ is *sequences on the same DNA
strands*, or in technical terms: sequences that are on the same feature,
regardless of strand orientation.

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
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Examples

``` r
gb       <- GRanges(c("Ref:100-200:+",   "Ref:400-500:+",    "Ref:600-700:+"))
gb$query <- GRanges(c("Que:1100-1200:+", "Que2:1700-1800:+", "Que:1500-1600:+"))
correlation_index(gb)
#> [1] 1
```
