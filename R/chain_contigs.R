#' Chain contigs in pairwise alignments
#'
#' Description TBW
#'
#' @param gb [`GBreaks`] object of the pairwise alignment.
#' @param tol Unaligned region of width lesser than or equal to `tol` in both
#'        the reference and query case will be bridged in coalescing.
#'
#' @return TBW
#'
#' @family Colinearity functions
#'
#' @author Charles Plessy
#'
#' @examples
#'
#' gb4       <- GenomicRanges::GRanges(c("Ref:100-200:+", "Ref:400-500:+", "Ref:600-700:+"))
#' gb4$query <- GenomicRanges::GRanges(c("Que:1100-1200:+", "Que:1700-1800:+", "Que:1500-1600:+"))
#' chain_contigs(gb4)
#'
#' @importFrom BiocParallel bplapply
#' @export

chain_contigs <- function(gb, tol = Inf, drop.within = FALSE) {
  # The rest of the algorithm assumes that the reference ranges are sorted
  gb <- sort(gb, ignore.strand = TRUE)

  # Split per strand, target name and query name
  gbl <- split(gb, paste(strand(gb), seqnames(gb), seqnames(gb$query)), drop = TRUE)
  chains <- BiocParallel::bplapply(gbl, coalesce_contigs, tol = tol) |>
    GRangesList() |> unlist() |> sort(ignore.strand = TRUE)

  if(!drop.within) return(chains |> unname())

  # Discard the chains that are included in other chains
  ovl <- findOverlaps(chains, type = "within")
  toRemove <- queryHits(ovl[width(chains[queryHits(ovl)]) < width(chains[subjectHits(ovl)])])
  chains[-queryHits(ovl)[toRemove]] |> unname()
}
