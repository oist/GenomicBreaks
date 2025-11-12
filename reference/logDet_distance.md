# LogDet (log-determinant) distance

Computes the LogDet (a.k.a. log-determinant) nucleotide distance between
two sequences from a 4×4 A/C/G/T joint count table. The method
normalizes the joint frequency matrix by the marginal base compositions
of the two sequences and then applies the log-determinant transform.

## Usage

``` r
logDet_distance(m, pseudocount = 0)
```

## Arguments

- m:

  A numeric matrix of pairwise counts with row/column names including
  `A,C,G,T`. If `m` is a list, the field `$probability_matrix` is used
  (compatible with your other helpers).

- pseudocount:

  A non-negative number added to each A/C/G/T cell before normalization
  (default `0`). Use a small value (e.g. `0.5`) to reduce the chance of
  a singular matrix when some joint cells are zero.
  [2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)

## Value

A single numeric: the LogDet distance (substitutions/site, unitless).
Returns `Inf` if any required determinant/sign is non‑positive.

A single numeric: the HKY85 distance (substitutions/site).

## Details

Let \\M\\ be the 4×4 **joint count** matrix over `A,C,G,T` and \\F = M /
\sum M\\ the **joint frequency** (divergence) matrix. Let
\\\boldsymbol{r} = F \mathbf{1}\\ (row sums, seq1 base frequencies) and
\\\boldsymbol{c} = F^\top \mathbf{1}\\ (column sums, seq2 base
frequencies). Define diagonal matrices \\D_R =
\mathrm{diag}(\boldsymbol{r})\\ and \\D_C =
\mathrm{diag}(\boldsymbol{c})\\. The **LogDet** normalization uses \$\$
S \\=\\ D_R^{-1/2}\\ F \\ D_C^{-1/2}, \$\$ and the **distance** is \$\$
d\_{\mathrm{LogDet}} \\=\\ -\\\ln \\\big( \det S \big) \\=\\ -\\\ln
\\\big( \det F \big) \\+\\ \tfrac{1}{2}\\\big\[ \ln \\\big( \det D_R
\big) \\+\\ \ln \\\big( \det D_C \big) \big\]. \$\$ This formulation
follows the LogDet/paralinear family where determinants commute (real
scalars), enabling robustness to composition heterogeneity; see Lockhart
et al. (1994) for the transform and discussion, and MEGA’s note on
non-computability when log terms approach zero.
[1](https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf)[2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)

## References

Lockhart PJ, Steel MA, Hendy MD, Penny D (1994). Recovering evolutionary
trees under a more realistic model of sequence evolution. *Mol Biol
Evol* **11**:605–612.
<https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf>
[1](https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf)

MEGA help: “LogDet Distance Could Not Be Computed” (notes on log terms).
<https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm>
[2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)

Gu X, Li WH (1996). Bias-corrected paralinear and LogDet distances and
tests under nonstationary nucleotide frequencies. *Mol Biol Evol*
**13**:1375–1383.
<https://paperity.org/p/54207439/bias-corrected-paralinear-and-logdet-distances-and-tests-of-molecular-clocks-and>
[3](https://paperity.org/p/54207439/bias-corrected-paralinear-and-logdet-distances-and-tests-of-molecular-clocks-and)

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
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`exampleSubstitutionMatrix`](https://oist.github.io/GenomicBreaks/reference/exampleSubstitutionMatrix.md),
[`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md)

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
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Examples

``` r
logDet_distance(exampleSubstitutionMatrix)
#> [1] 1.130934
```
