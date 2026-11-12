# Monkey indices

Twelve summary statistics of a pairwise genome alignment that were not
designed to represent a biological fact. They were retained from a
screen of 5,420 functions assembled mechanically from the columns of the
alignment table and from generic descriptive statistics ("monkey
typing"), for the following properties on 1,910 alignments of 16 clades:
reproducible between the alignment of *A* to *B* and the alignment of
*B* to *A*, more alike between genome pairs that share an ancestor than
between random pairs, different from the value obtained after shuffling
the order and the strand of the alignments, and not strongly correlated
with the statistics that already existed in this package, nor with each
other (absolute rank correlation below 0.6). Their biological meaning,
if any, is to be found. Each function is named after a famous monkey.

## Usage

``` r
wukong_index(gb)

abu_index(gb)

albert_index(gb)

marcel_index(gb)

imo_index(gb)

rafiki_index(gb)

hanuman_index(gb)

jack_index(gb)

able_index(gb)

mizaru_index(gb)

george_index(gb)

baker_index(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object, not coalesced, in which all ranges are on the `+` or `-`
  strand.

## Value

A numeric value, or `NA` when it can not be computed: fewer than three
values, or no variation.

## Details

All the functions take the alignments as loaded and sort them on the
*target* genome. `mizaru_index()` then coalesces them with
[`coalesce_contigs`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md)
and its default arguments; the others use them as they are. They use the
following quantities.

- A *junction* is a pair of consecutive alignments on the same *target*
  sequence. A *same-query junction* is a junction whose two alignments
  are also on the same *query* sequence.

- The *target gap* of a junction is the number of bases between its two
  alignments on the *target* genome, plus one; the *query gap* of a
  same-query junction is the same on the *query* genome, whatever the
  order and the strand of the two alignments.

- The *gap contrast* of a same-query junction is \\(t - q) / (t + q)\\
  where \\t\\ and \\q\\ are its target gap and its query gap. It is
  positive when the two alignments are further apart on the *target*.

- The *rank step* of a same-query junction is the number of positions
  between its two alignments in the order of the alignments on the
  *query* genome (1 for neighbours), and its *direction* is +1 when the
  two alignments are on the same strand and follow each other on the
  *query* as that strand implies, -1 when they follow each other in the
  opposite direction, and 0 when they are on opposite strands.

- The *edge distance* of an alignment on a genome is the distance from
  its centre to the nearest end of the sequence, divided by half the
  length of the sequence: 0 at the ends, 1 in the middle. The *relative
  position* is the position of the centre from 0 (first base) to 1 (last
  base).

- *Sizes* (widths, gaps, rank steps) are used on the logarithmic scale.

The functions are the following. The monkey is introduced first.

- `wukong_index()`:

  *Sun Wukong: the Monkey King of the Chinese novel Journey to the
  West.* Tail weight of the query gaps: on the logarithmic scale, range
  between the 10th and the 90th percentiles divided by the interquartile
  range.

- `abu_index()`:

  *Abu: kleptomaniac pet monkey of the street thief in the Aladdin
  films.* Quartile skewness of the gap contrasts: \\(Q_3 + Q_1 - 2 Q_2)
  / (Q_3 - Q_1)\\.

- `albert_index()`:

  *Albert II: rhesus macaque, first primate to reach space, in 1949.*
  Excess kurtosis of the gap contrasts.

- `marcel_index()`:

  *Marcel: white-headed capuchin kept by Ross in the sitcom Friends.*
  Standard deviation of the absolute difference between the relative
  positions of an alignment on the two genomes. The *query* position is
  counted from the other end of the sequence when most aligned bases
  between the two sequences are on the minus strand.

- `imo_index()`:

  *Imo: Japanese macaque of Koshima who started washing sweet potatoes.*
  Correlation between the edge distance of an alignment on the *target*
  genome and the edge distance of its neighbours on the *query* genome.
  Neighbours are the alignments before and after it on the same *target*
  sequence; the two members of a junction are exchangeable.

- `rafiki_index()`:

  *Rafiki: mandrill shaman who presents the newborn cub in The Lion
  King.* Rank correlation, over same-query junctions, between the query
  gap and the mean edge distance of the two alignments on the *query*
  genome.

- `hanuman_index()`:

  *Hanuman: monkey-faced Hindu deity, devoted companion of Rama in the
  Ramayana.* Rank correlation, over junctions, between the mean width of
  the two alignments and the absolute difference of their edge distances
  on the *query* genome.

- `jack_index()`:

  *Jack: chacma baboon who operated railway signals in South Africa in
  the 1880s.* Rank correlation, over same-query junctions, between the
  rank step and the direction.

- `able_index()`:

  *Able: rhesus macaque, Miss Baker's companion on the 1959
  spaceflight.* Correlation, over same-query junctions, between the rank
  step and the absolute difference of the edge distances of the two
  alignments on the *query* genome.

- `mizaru_index()`:

  *Mizaru: the see-no-evil monkey carved at the Toshogu shrine in
  Nikko.* Rank correlation, over the junctions of the coalesced
  alignments, between the target gap and the absolute difference of the
  edge distances of the two coalesced alignments on the *query* genome.

- `george_index()`:

  *Curious George: inquisitive little monkey of the children's books by
  the Reys.* Correlation, over junctions, between the target gap and the
  mean score per aligned column of the two alignments. Stops with an
  error without the `score` and `aLength` metadata columns.

- `baker_index()`:

  *Miss Baker: squirrel monkey who survived a suborbital spaceflight in
  1959.* Correlation, weighted by the width of the alignments on the
  *target* genome, between the number of gap openings per aligned column
  of an alignment and the geometric mean of the target gaps of the
  junctions that flank it. Stops with an error without the `gapOpens`
  and `aLength` metadata columns.

Ranks are computed within each vector before restricting it to the
junctions used: in `rafiki_index()`, the mean edge distances are ranked
among all junctions. Variances are population variances. Quantiles are
those of [`stats::quantile`](https://rdrr.io/r/stats/quantile.html) with
default settings.

The results do not depend on the order of the sequences nor on their
orientation: reversing a sequence with
[`reverse`](https://oist.github.io/GenomicBreaks/reference/reverse.md)
does not change them. They are not symmetric: the value on `swap(gb)`
differs. Note that
[`swap`](https://oist.github.io/GenomicBreaks/reference/swap.md) drops
the metadata columns that `george_index()` and `baker_index()` need.

The functions that use edge distances or relative positions need
sequence lengths. Missing lengths are replaced by the end of the last
alignment of each sequence, with
[`forceSeqLengths`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md).

## See also

The vignette
[`vignette("MonkeyIndices", package = "GenomicBreaks")`](https://oist.github.io/GenomicBreaks/articles/MonkeyIndices.md)
reports how the functions were selected and how they behave on the
alignments of the screen. The code of the screen is in
`system.file("scripts", "monkey_screen", package = "GenomicBreaks")`.

## Author

Claude (Anthropic), from a screen requested by Charles Plessy.

## Examples

``` r
# Alignment of two yeast genomes
gb <- load_genomic_breaks(
  system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"))
wukong_index(gb)
#> [1] 2.496076
abu_index(gb)
#> [1] 0.1227582
albert_index(gb)
#> [1] -1.082118
marcel_index(gb)
#> [1] 0.1281138
imo_index(gb)
#> [1] 0.9163483
rafiki_index(gb)
#> [1] 0.03985051
hanuman_index(gb)
#> [1] 0.5777603
jack_index(gb)
#> [1] -0.5166595
able_index(gb)
#> [1] 0.2710925
mizaru_index(gb)
#> [1] 0.007696371

# The orientation of the sequences does not matter.
marcel_index(reverse(forceSeqLengths(gb)))
#> [1] 0.1281138

# The GFF3 format does not record the alignment length nor the number of gap
# openings, which the blasttab+ format does.  Lower bounds are used here
# for illustration.
gb$aLength  <- pmax(width(gb), width(gb$query))
gb$gapOpens <- as.integer(width(gb) != width(gb$query))
george_index(gb)
#> [1] 0.1178123
baker_index(gb)
#> [1] -0.01414999
```
