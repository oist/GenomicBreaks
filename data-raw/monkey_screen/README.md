# Screen that produced the monkey indices

These scripts screened 5,420 mechanically assembled summary functions of a
pairwise genome alignment and retained the ten that are exported by the
package as `?monkey_indices`.  They are kept for the record: they do not
regenerate anything in `data/` or `inst/extdata/`.

The screen was run on 2026-11-11 on 1,910 alignments (955 genome pairs in both
orders, 182 genomes, 16 clades) of the project _Scrambling in the Tree of
Life_.  Python 3.13 with numpy, scipy, pandas, scikit-learn and PyYAML; two
cores, 7 GB of memory, about 2.5 hours of computation.

## Input

`extractMonkeySubset.py` selects the genomes and writes one tar file per
Branch.  For each Branch that has both alignments
(`AlignedGenomes/3.0/<seed>/<Branch>`) and statistics
(`StatsOnPairs/latest/<Branch>`), it computes the P distance of every ordered
pair from the `.o2o.matrix.txt` files, builds a UPGMA tree, cuts it in
`--clades` clades (2 in the run that was used), picks `--genomes` genomes (12)
and exports, for every ordered pair of them, the `.blasttabplus.gz`,
`.o2o.matrix.txt`, `.o2o.tsv` and `.yaml` files.  The tar files are unpacked
in `$MONKEY_ROOT/data/monkey/`.

## Steps

Set `MONKEY_ROOT` to a working directory that contains `data/monkey/`, and run
the scripts in order.  `mk_*.py` are modules.

| Script | Does | Time |
|---|---|---|
| `01_build_cache.py` | parses the alignments once, applies the filter of the statistics pipeline (250 sequences with the most aligned bases) | 3 min |
| `02_existing.py` | table of the YAML statistics | 1 min |
| `03_validate.py` | checks the loader and the coalescing against the YAML counts and totals | 1 min |
| `04_test_engine.py` | checks the candidate engine on four alignments | 1 min |
| `05_stage1.py` | all candidates on all alignments, as they are and with target and query exchanged | 50 min |
| `06_extras.py` | adds to the reference table the package indices that are not in the YAML files, and simple derivations | 3 min |
| `07_check_existing.R` | values of the R package, to compare with `mk_existing.py` | 1 min |
| `08_metrics.py` | definedness, reciprocal agreement, ancestor-node agreement, correlation with the references | 5 min |
| `09_invariance.py` | same value when sequences are reversed and renamed | 1 min |
| `10_null.py` | all candidates on 5 shuffles of each genome pair | 65 min |
| `11_redundancy.py` | predictability from all references together (ridge, leave one Branch out) | 2 min |
| `12_shortlist.py` | joins the metrics (`results/metrics_all.tsv`) | 1 min |
| `13_null_metrics.py` | departure from the shuffles | 1 min |
| `14_select.py` | filters, relaxation, greedy selection of 10 (`results/selected.tsv`) | 2 min |
| `15_final.py` | recomputes the selected candidates on every alignment (`results/final_values.tsv`) | 5 min |
| `16_r_check.R` | values of the R functions on a list of alignments, also after `swap()` and `reverse()` | 2 min |
| `17_final_metrics.py` | all metrics of the selected candidates, from the final values | 2 min |

## Candidates

A candidate is named `<level>.<reducer>.<vector>[~<vector>]`: level `raw`
(alignments as loaded) or `coa` (coalesced); see the header of `mk_engine.py`
for the vectors and the reducers.  The exported functions are:

| R function | Candidate |
|---|---|
| `wukong_index()` | `raw.tailw.qstep` |
| `abu_index()` | `raw.bowley.gdiff` |
| `albert_index()` | `raw.kurt.gdiff` |
| `marcel_index()` | `raw.sd.pdisp` |
| `imo_index()` | `raw.xcor.tedge~qedge` |
| `rafiki_index()` | `raw.rcor.qstep~jm_qedge` |
| `hanuman_index()` | `raw.rcor.jm_tw~jd_qedge` |
| `jack_index()` | `raw.rcor.rstep~dirc` |
| `george_index()` | `raw.pcor.tgap~jm_sdens` |
| `baker_index()` | `raw.wcor.gapf~bgap` |

## Filters

The docstring of `14_select.py` lists the filters and the two settings that
were relaxed until ten functions remained: the common ceiling on correlations
(0.6 in the end) and the share of genome pairs in which a function must depart
from the shuffles (one half in the end).

## Notes

- The stage-1 values were computed before late changes of `mk_engine.py`
  (ties of the majority strand; a zero-variance guard).  The selected
  candidates were recomputed afterwards: one of them, `raw.sd.pdisp`, changed
  by at most 0.13 %.
- The shuffles of the coalesced level shuffle the coalesced blocks and do not
  coalesce them again.
- No alignment of this data set has more than 250 sequences, so the filter of
  the statistics pipeline was never applied.
