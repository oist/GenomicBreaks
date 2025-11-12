# Hasegawa–Kishino–Yano (1985) distance (HKY85)

Computes the HKY85 nucleotide substitution distance between two
sequences, using only the A/C/G/T block of a pairwise count matrix. Base
frequencies are estimated as the average of row and column marginals,
and transitions are pooled into a single class (A↔G and C↔T).

## Usage

``` r
HKY85_distance(m)
```

## Arguments

- m:

  A numeric matrix of pair counts with row/column names including
  `A,C,G,T`. If `m` is a list, the field `$probability_matrix` is used
  (for compatibility with your other helpers).

## Value

A single numeric: the HKY85 distance (substitutions/site).

## Details

Let \\P\\ be the total fraction of **transition** differences (A↔G
**plus** C↔T) among A/C/G/T pairs, and \\Q\\ the fraction of
**transversion** differences (all other A↔{C,T}, G↔{C,T} pairs). Let
\\\pi_A,\pi_C,\pi_G,\pi_T\\ be the **average base frequencies** across
the two sequences (row and column marginals averaged), \\\pi_R = \pi_A +
\pi_G\\, and \\\pi_Y = \pi_C + \pi_T\\.

The HKY85 distance (special case of TN93 with a single transition rate)
is: \$\$ d\_{\mathrm{HKY85}} = -\\2(\pi_A \pi_G + \pi_C \pi_T)\\
\ln\\\left(1 - \frac{P}{2(\pi_A \pi_G + \pi_C \pi_T)}\right) -\\2\\\pi_R
\pi_Y\\ \ln\\\left(1 - \frac{Q}{2\\\pi_R \pi_Y}\right). \$\$

Implementation notes:

- Rows/columns outside `A,C,G,T` are dropped.

- \\P\\ and \\Q\\ are computed from the normalized 4×4 table.

- \\\pi\_\bullet\\ are estimated as the average of row and column
  marginals.

- Returns `Inf` if any log argument is non‑positive or a denominator is
  zero.

## References

Hasegawa, M., Kishino, H., and Yano, T. (1985). Dating of the human–ape
splitting by a molecular clock of mitochondrial DNA. *J. Mol. Evol.*
22:160–174.

RevBayes tutorial (overview of HKY, parameters \\\kappa,\pi\\):
<https://revbayes.github.io/tutorials/ctmc/>

IQ‑TREE model manual (relationships among HKY, K80, TN93):
<https://iqtree.github.io/doc/Substitution-Models>

## See also

Other Alignment statistics:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
[`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md),
[`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
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
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Entirely written by Copilot and not proofchecked yet

## Examples

``` r
HKY85_distance(exampleSubstitutionMatrix)
#> [1] 0.288648
```
