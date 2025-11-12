# Index measuring how 'random' the alignment strand has become.

In groups of species where major changes of gene order happened but
genes tended to stay on the same chromosome, a marked feature is that
the strand on which homologous sequences align appears to be randomised.
This index expresses it with a number.

## Usage

``` r
strand_randomisation_index(gb, tiles = NULL)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- tiles:

  A number of tiles

## Value

Returns a numeric value between 0 and 1.

## Details

For each feature in the *target* genome, the total length of alignments
on the plus strand is subtracted from that on the minus strand, and the
absolute value is taken (because there is no guarantee that homologous
chromosomes are sequenced in the same direction in different
assemblies). This value is normalized by the total number of aligned
bases. A weighted mean is then computed across all features, using each
feature’s total aligned bases as its weight. Thus, a number close to 1
is expected for closely related genomes.

## Note

This index is designed for comparison of chromosomal assemblies that
have a strong conservation of synteny in the sense most homologous genes
are mapped on pairs of homologous chromosomes. In 2025, a new `tiles`
option was added to make the index more robust to assembly errors and
karyotype changes.

## References

Plessy C, Mansfield MJ, Bliznina A, Masunaga A, West C, Tan Y, Liu AW,
Grašič J, Del Río Pisula MS, Sánchez-Serna G, Fabrega-Torrus M,
Ferrández-Roldán A, Roncalli V, Navratilova P, Thompson EM, Onuma T,
Nishida H, Cañestro C, Luscombe NM (2024). Extreme genome scrambling in
marine planktonic *Oikopleura dioica* cryptic species. *Genome
Research*, 34(3), 426–440. <https://doi.org/10.1101/gr.278295.123>.
PMID: 38621828

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
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Charles Plessy

## Examples

``` r
strand_randomisation_index(exampleInversion)
#> [1] 0.3992674
```
