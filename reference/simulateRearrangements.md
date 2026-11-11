# Simulate rearranged genomes

Simulates a pair of genomes made of `n` aligned blocks: the *target*
genome keeps the blocks in their ancestral order, and the *query* genome
is derived from it by random inversions and transpositions. Both genomes
are then cut into sequences at random positions.

## Usage

``` r
simulateRearrangements(
  n = 14,
  nEvents = 3,
  types = c("inversion", "transposition", "inverted transposition"),
  maxMove = 2,
  nTargetSeqs = 1,
  nQuerySeqs = 1,
  ...
)
```

## Arguments

- n:

  Number of blocks.

- nEvents:

  Number of events.

- types:

  Types of events to sample from (see Details).

- maxMove:

  Maximum number of blocks moved by a transposition.

- nTargetSeqs, nQuerySeqs:

  Number of sequences in each genome.

- ...:

  Passed to
  [`permutationVectorToGBreaks()`](https://oist.github.io/GenomicBreaks/reference/permutationVectorToGBreaks.md)
  (`width`, `spacing`).

## Value

A sorted
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object of length `n`. Its `metadata()` list contains an `events` data
frame with one row per event, in order: `type`; `segment`, the blocks
affected by the event, in *query* order before the event, with a minus
sign for the blocks on the minus strand; and `after`, for
transpositions, the block preceding the moved segment in the *query*
genome after the event (`0` when the segment was moved to the start).
Blocks are identified by their position in the *target* genome, which is
the ancestral order.

## Details

Each event picks its blocks at random positions of the *query* genome as
it stands after the previous events, before it is cut into sequences:

- `inversion`: reverses a segment of one block or more and flips its
  strand.

- `transposition`: moves a segment of `1` to `maxMove` blocks to another
  position.

- `inverted transposition`: same, and inverts the segment.

Events are not checked for parsimony: a later event can undo, or extend,
the effect of an earlier one. The random cuts of each genome into
sequences are not recorded as events.

## See also

[`permutationVectorToGBreaks()`](https://oist.github.io/GenomicBreaks/reference/permutationVectorToGBreaks.md)
and
[`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md).
A set of 200 simulations is provided in the
`simulatedRearrangements.rds` file of the package (see Examples).

## Examples

``` r
set.seed(1)
gb <- simulateRearrangements(n = 10, nEvents = 2)
gb
#> GBreaks object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |           query
#>           <Rle> <IRanges>  <Rle> |       <GRanges>
#>    [1]    chrT1   100-189      + |   chrQ1:100-189
#>    [2]    chrT1   200-289      - |   chrQ1:700-789
#>    [3]    chrT1   300-389      - |   chrQ1:600-689
#>    [4]    chrT1   400-489      + |   chrQ1:200-289
#>    [5]    chrT1   500-589      + |   chrQ1:300-389
#>    [6]    chrT1   600-689      + |   chrQ1:400-489
#>    [7]    chrT1   700-789      + |   chrQ1:500-589
#>    [8]    chrT1   800-889      + |   chrQ1:800-889
#>    [9]    chrT1   900-989      + |   chrQ1:900-989
#>   [10]    chrT1 1000-1089      + | chrQ1:1000-1089
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
metadata(gb)$events
#>        type         segment after
#> 1 inversion         4 5 6 7    NA
#> 2 inversion 2 3 -7 -6 -5 -4    NA
permutationVector(gb)
#>  [1]  1 -7 -6  2  3  4  5  8  9 10
flagAll(gb)
#> GBreaks object with 10 ranges and 2 metadata columns:
#>        seqnames    ranges strand |           query        flag
#>           <Rle> <IRanges>  <Rle> |       <GRanges> <character>
#>    [1]    chrT1   100-189      + |   chrQ1:100-189        <NA>
#>    [2]    chrT1   200-289      - |   chrQ1:700-789         Col
#>    [3]    chrT1   300-389      - |   chrQ1:600-689        <NA>
#>    [4]    chrT1   400-489      + |   chrQ1:200-289         Col
#>    [5]    chrT1   500-589      + |   chrQ1:300-389         Col
#>    [6]    chrT1   600-689      + |   chrQ1:400-489         Col
#>    [7]    chrT1   700-789      + |   chrQ1:500-589        <NA>
#>    [8]    chrT1   800-889      + |   chrQ1:800-889         Col
#>    [9]    chrT1   900-989      + |   chrQ1:900-989         Col
#>   [10]    chrT1 1000-1089      + | chrQ1:1000-1089        <NA>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Pre-computed simulations, stored as permutation vectors
sims <- readRDS(system.file("extdata/simulatedRearrangements.rds",
                            package = "GenomicBreaks"))
sims[[1]]
#> $v
#>  [1]  1  2  3  5  6  7  8  9 10 -4 11 12 13 14
#> 
#> $tsplit
#> [1]  4 10
#> 
#> $qsplit
#> [1] 14
#> 
#> $events
#>                     type segment after
#> 1 inverted transposition      10     3
#> 
#> $params
#> $params$n
#> [1] 14
#> 
#> $params$nEvents
#> [1] 1
#> 
#> $params$nTargetSeqs
#> [1] 2
#> 
#> $params$nQuerySeqs
#> [1] 1
#> 
#> 
with(sims[[1]], permutationVectorToGBreaks(v, tsplit, qsplit))
#> GBreaks object with 14 ranges and 1 metadata column:
#>        seqnames    ranges strand |           query
#>           <Rle> <IRanges>  <Rle> |       <GRanges>
#>    [1]    chrT1   100-189      + |   chrQ1:100-189
#>    [2]    chrT1   200-289      + |   chrQ1:200-289
#>    [3]    chrT1   300-389      + |   chrQ1:300-389
#>    [4]    chrT1   400-489      + |   chrQ1:500-589
#>    [5]    chrT2   100-189      + |   chrQ1:600-689
#>    ...      ...       ...    ... .             ...
#>   [10]    chrT2   600-689      - |   chrQ1:400-489
#>   [11]    chrT2   700-789      + | chrQ1:1100-1189
#>   [12]    chrT2   800-889      + | chrQ1:1200-1289
#>   [13]    chrT2   900-989      + | chrQ1:1300-1389
#>   [14]    chrT2 1000-1089      + | chrQ1:1400-1489
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
