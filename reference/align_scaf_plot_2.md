# Compare Alignments over two Reference scaffolds

The function plots the pairwise alignments from two of the reference
scaffolds. This can be used as an aid in super-scaffolding. X-axis is
the position on the reference scaffolds, and y-axis is aligned query
scaffold.

## Usage

``` r
align_scaf_plot_2(gr_ob, scaf)
```

## Arguments

- gr_ob:

  GRanges object containing pairwise alignment

- scaf:

  A list containing two character objetcs of the names of the reference
  scaffolds you wish to compare on the plot

## Value

scaffold alignment plot, comparing two reference scaffolds
