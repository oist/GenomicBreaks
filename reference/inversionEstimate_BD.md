# Inversion Estimate - Berestycki and Durrett (2006)

Computes the expected number of inversions that explains the differences
between the gene orders of two genomes.

## Usage

``` r
inversionEstimate_BD(gb, mode = "few", chrom_sep = "---")
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- mode:

  Optional parameter that specifies how estimates are computed. Two
  options are available:

  1.  `few` : recommended when computing a few values in one function
      call;

  2.  `many` : use when computing many values.

- chrom_sep:

  An optional string that separates the chromosome names in each item of
  the returned list.

## Value

A list where each item corresponds to a chromosome pair. The name of
each item is the matched reference and query chromosome names
concatenated, separated by `chrom_sep`. Each item contains the following
information: 1. Basic stats of the breakpoint graph: `N`,
`nbBreakpoints`, `nbCycles`. See the documentation of
`breakpointGraphProperties` for more details; 2. Details of the
*Berestycki and Durrett* estimate stored in `expinv_BD`: \* `k_beg` and
`k_end`: the lower and upper bounds for the expected number of
inversions, respectively; \* `k_avg` : the average of `k_beg` and
`k_end`; \* `nb_cyc_beg` and `nb_cyc_end`: the expected number of cycles
in the breakpoint graph after `k_beg` and `k_end` inversions. The
observed number of cycles in the breakpoint graph of the two genomes
should fall between those values.

## Details

To compute the expected number of inversions, the method needs to first
compute the *breakpoint graph*, a classical data structure introduced by
Hannenhalli and Pevzner and used in many genome rearrangement problems.

Given the number of cycles in the breakpoint graph, it returns the
estimated number of evolutionary events (inversions).

This function implements Theorem 3 from Berestycki and Durrett (2006),
also described in Equation 6 in Biller et al. (2015).

## Note

If the input genomes are multichromosomal, it is recommended to call the
function
[`matchPairs`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md)
before this function, to create a mapping between reference and query
chromosomes.

## References

Berestycki, Nathanaël, and Rick Durrett. "A phase transition in the
random transposition random walk." Discrete Mathematics and Theoretical
Computer Science. Discrete Mathematics and Theoretical Computer Science,
2003.

Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into
turnip: polynomial algorithm for sorting signed permutations by
reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.

Biller, Priscila, Laurent Guéguen, and Eric Tannier. "Moments of genome
evolution by double cut-and-join." BMC bioinformatics 16.Suppl 14
(2015): S7.

## See also

[`breakpointGraphProperties`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md)
for computing key properties of the breakpoint graph needed for this
estimate.

Other Rearrangement distances:
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`JC69_distance_allseq()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance_allseq.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md),
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
if (FALSE) { # \dontrun{
# Create a chromosome mapping given a GBreaks object (useful when genomes are multichromosomal).
chrMapping <- matchPairs(exampleInversionBader2001)
# Compute the expected number of inversions using the method from Berestycki and Durrett (2006).
expNbInversions <- inversionEstimate_BD(chrMapping)
# Compute the minimum number of inversions using the method from Hannehalli and Pevzner (1999).
minNbInversions <- inversionDistance(chrMapping)
# Output: Example from Bader et al. (2001): 
#         - Inversion distance = 7 , Expected nb. of inversions = 7
cat("Example from Bader et al. (2001):")
cat(paste("- Inversion distance =", minNbInversions))
cat(paste("- Expected nb. of inversions =", expNbInversions[[1]]$expinv_BD$k_avg))

# Another example, this time without computing the chromosome mapping.
# The chromosome mapping is not needed if genomes are unichromosomal.
expNbInversions <- inversionEstimate_BD(exampleInversionBergeron2005b)
minNbInversions <- inversionDistance(exampleInversionBergeron2005b)
# Output: Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6): 
#         - Inversion distance = 13 , Expected nb. of inversions = 15
cat("Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6):")
cat(paste("- Inversion distance =", minNbInversions))
cat(paste("- Expected nb. of inversions =", expNbInversions[[1]]$expinv_BD$k_avg))
} # }
```
