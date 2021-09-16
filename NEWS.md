# GenomicBreaks 0.7.0

* New `mergeSeqLevels()` function to concatenate sequence features, for instance
  before plotting.
* New `reverse()` method for `GBreaks` and `GRanges` objects.
* Really return a `GBreaks` object in `coalesce_contigs()`.
* Pass extra arguments from `plotApairOfChrs()` to `genoPlotR::plot_gene_map`.

# GenomicBreaks 0.6.0

* Add pairwise chromosome plot functions using the
  [`genoPlotR`](http://genoplotr.r-forge.r-project.org) package,
  see in particular `plotApairOfChrs()`.

# GenomicBreaks 0.5.3

* _pkgdown_ update.

# GenomicBreaks 0.5.2

* Greatly speed up `flagLongShort()`.

# GenomicBreaks 0.5.1

* Correct a variable name bug in `flagLongShort()`.
* Export the `GBreaks()` constructor.

# GenomicBreaks 0.5.0

* Add a `flagLongShort()` function to transfer annotations about chromosome arms.

# GenomicBreaks 0.4.0

* Use alignments of _Saccharomyces paradoxus_ to _S. cerevisiae_ as example data.
* Addition of a `GOC()` function to calculate Gene Order Conservation.

# GenomicBreaks 0.3.0

* Added the `leftInversionGaps()` function.

# GenomicBreaks 0.2.0

* Document the package with [`pkgdown`](https://pkgdown.r-lib.org/).
