suppressMessages(library(GenomicBreaks))
args <- commandArgs(trailingOnly = TRUE)
for (f in args) {
  gb <- load_genomic_breaks(f)
  sw <- sort(swap(gb), ignore.strand = TRUE)
  bp <- function(gb) {
    l <- split(gb, paste(seqnames(gb), seqnames(gb$query), sep = "---"), drop = TRUE)
    r <- lapply(l, breakpointGraphProperties)
    c(sum(sapply(r, `[[`, "nbBreakpoints")), sum(sapply(r, `[[`, "N")) - sum(sapply(r, `[[`, "nb_cycles")))
  }
  cat(basename(f), length(gb), karyotype_index(gb), karyotype_index(sw), tau_index(gb), tau_index(sw),
      bp(gb), synteny_index(gb), correlation_index(gb), strand_randomisation_index(gb), GOC(gb), "\n")
}
