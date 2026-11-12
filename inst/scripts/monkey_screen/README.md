# Screen that produced the monkey indices

Python code of the screen that selected the twelve functions exported by the
package as `?monkey_indices`.  It is kept for the record: no function of the
package uses it, and it does not generate anything in `data/` or
`inst/extdata/`.

The vignette _Monkey indices_ (`vignette("MonkeyIndices", package =
"GenomicBreaks")`, source in `vignettes/MonkeyIndices.Rmd`) describes the
alignments, the candidate functions, the filters and the results, and lists
what each script does.

To run the screen: Python 3 with numpy, scipy, pandas, scikit-learn and
PyYAML.  Set the environment variable `MONKEY_ROOT` to a working directory
that contains the alignments in `data/monkey/` (the tar files written by
`extractMonkeySubset.py`, unpacked), and run the numbered scripts in order.
The `mk_*.py` files are modules.

`19_vignette_tables.py <directory>` writes the two tables that the vignette
reads, `MonkeyIndices_values.tsv` and `MonkeyIndices_functions.tsv`, which are
kept in `vignettes/`.
