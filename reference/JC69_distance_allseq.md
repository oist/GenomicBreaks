# Variation of the Jukes-Cantor 1969 distance

The Jukes-Cantor 1969 (JC69) distance corrects the p-distance for
multiple substitutions, providing an estimate of evolutionary distance
that is proportional to time under the model. The correction is based on
the proportion of nucleotide differences, typically obtained by counting
mismatches in the alignment matrix. Here, the Jukes-Cantor equation
remains unchanged, but the definition of which base pairs are considered
different is modified.

## Usage

``` r
JC69_distance_allseq(gb, m, adjust_p = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- m:

  A matrix of **counts** for bases of the *target* genome to be aligned
  to bases on the *query* genome.

- adjust_p:

  A boolean flag. If `TRUE`, the distance is scaled between `0` and
  `0.75` to ensure the logarithm stays positive.

## Value

Returns a numeric value representing the evolutionary distance between
two genomes. The greater the value, the more genetically different the
genomes are.

## Details

In this function, the fraction of nucleotides that are different
incorporates not only the mismatches from the alignment matrix but also
the base pairs that are left unaligned. The rationale is that these
unaligned base pairs likely differ primarily because of point
substitutions and should therefore be treated as mismatches that were
not detected by the aligner. In fact, the aligner can only spot
mismatches in regions where the proportion of matches is high enough for
alignment.

Notice that gaps are usually not included in the Jukes-Cantor distance,
including in this variation, because they are generally considered to
result from indels (insertions and deletions that affect multiple
nucleotides at once), whereas the Jukes-Cantor model is based only on
point substitutions. Therefore, including regions that were likely
affected by large evolutionary events, such as gaps caused by indels,
would incorrectly inflate a distance calculated under a model in which
only one position is mutated at a time.

## References

Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules." In
*Mammalian Protein Metabolism* (pp. 21–132). Academic Press.

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
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
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
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

Priscila Biller

## Examples

``` r

# Only the sequence length is used from the GenomicBreaks object.
gb       <- GRanges(c("Ref:100-35000000:+"))
gb$query <- GRanges(c("Que:1100-35000500:+"))
d <- JC69_distance_allseq(gb, exampleSubstitutionMatrix)
```
