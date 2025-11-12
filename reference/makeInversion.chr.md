# Make an inversion

Make an inversion

## Usage

``` r
makeInversion.chr(gr, dist = runif)
```

## Arguments

- gr:

  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object

- dist:

  Distribution function for new breakpoints on `gr`'s sequence levels.

## Examples

``` r
gr <- GRanges("chr1:1-1000000:+")
gr2 <- GRanges(c("chr1:1-1000000:+", "chr2:1-1000000:+"))
gr3 <- GRanges(c("chr1:1-1000000:+", "chr1:1000000-2000000:+"))
gr4 <- GRanges(c("chr1:1-1000000:+", "chr1:1000000-2000000:-"))
gr5 <- GRanges(c("chr1:1-1000000:-"))
```
