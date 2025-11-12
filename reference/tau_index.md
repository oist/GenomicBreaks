# Tau index

Tissue-specificity index from Yanai and coll., 2005 repurposed in 2022
by Ranz and coll. for synteny analysis.

## Usage

``` r
tau_index(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns a numeric value between 0 and 1.

## Details

In 2022, Ranz and coll., the tau index was computed for each chromosome
in a *target* genome compared to a *query*, with a formula such as: for
a given *target* feature, count its one-to-one orthologues on every
feature \\i\\ of the *query* genome, normalize these counts by dividing
by the largest count, subtract each normalized value from one, sum the
results, and divide by the number \\n\\ of *query* features minus one.

\$\$\frac{1}{n-1}\sum\_{i=1}^{n}\left(1 - \frac{x_i}{\max(x)}\right)\$\$

Here, the index is extended to the whole genome comparisons by computing
it for each feature of the *target* genome and returning the average
weighted by feature length.

## Note

Note that calculating the tau index on whole-genome nucleotide
alignments is not expected to produce meaningful results. This function
is more useful when comparing the positions of protein orthologues.

## References

Yanai I, Benjamin H, Shmoish M, Chalifa-Caspi V, Shklar M, Ophir R,
Bar-Even A, Horn-Saban S, Safran M, Domany E, Lancet D, Shmueli O.
(2005). Genome-wide midrange transcription profiles reveal expression
level relationships in human tissue specification. *Bioinformatics*,
21(5):650-659.
[doi:10.1093/bioinformatics/bti042](https://doi.org/10.1093/bioinformatics/bti042)
. PMID: 15388519

Ranz JM, González PM, Su RN, Bedford SJ, Abreu-Goodger C, Markow T.
(2022). Multiscale analysis of the randomization limits of the
chromosomal gene organization between Lepidoptera and Diptera. *Proc
Biol Sci*, 289(1967):20212183.
[doi:10.1098/rspb.2021.2183](https://doi.org/10.1098/rspb.2021.2183) .
PMID: 35042416

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
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md)

## Author

Charles Plessy

## Examples

``` r
tau_index(exampleTranslocation)
#> [1] 0.5
GenomicBreaks:::.tau_index(c(0,8,0,0,0,2,0,2,0,0,0))
#> [1] 0.95
```
