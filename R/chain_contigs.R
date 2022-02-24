#' Chain contigs in pairwise alignments
#'
#' Description TBW
#'
#' @param gb [`GBreaks`] object of the pairwise alignment.
#' @param tol Unaligned region of width lesser than or equal to `tol` in both
#'        the reference and query case will be bridged in coalescing.
#' @param drop.within Discard the chains that are included in other chains.
#'
#' @return TBW
#'
#' @family Colinearity functions
#'
#' @author Charles Plessy
#'
#' @examples
#' chain_contigs(exampleInversion)
#'
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors subjectHits
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
