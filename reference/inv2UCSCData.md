# Represent inversions as UCSCData objects

Scan a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object for inversions and return UCSC track object in
[`rtracklayer::UCSCData`](https://rdrr.io/pkg/rtracklayer/man/UCSCData-class.html)
format for export in BED12 format with a command such as
`rtracklayer::export(x, "test.bed", "BED")`.

## Usage

``` r
inv2UCSCData(gb)
```

## Arguments

- gb:

  A `GBreaks` object.

## Value

Returns a `UCSCData` object.
