#' Match pairs of sequences
#'
#' For each sequence of the _target_ genome, assign a sequence of the _query_
#' genome, and discard all the ranges that are not from matched pairs.
#'
#' @note This is not a reciprocal best match function, because it assumes that
#' the genomes are close enough and that there were no major karyotype changes,
#' so that the pairing is a trivial operation.
#'
#' @param gb A [`GBreaks`] object.
#' @param drop Drop unused sequence levels.
#'
#' @family modifier functions
#'
#' @return A `GBreaks` object.
#'
#' @author Charles Plessy
#'
#' @export
#'
#' @examples
#' matchPairs(exampleTranslocation)

matchPairs <- function(gb, drop = FALSE) {
  if (length(gb) == 0) return(gb)
  gbl <- split(gb, seqnames(gb))
  gb <- endoapply(gbl, \(gb){
    bestMatch <- tapply(width(gb), seqnames(gb$query), sum) |> sort() |> tail(1) |> names()
    gb[seqnames(gb$query) == bestMatch]
  }) |> unlist()
  if (isTRUE(drop)) {
    seqlevels(gb      ) <- seqlevelsInUse(gb      )
    seqlevels(gb$query) <- seqlevelsInUse(gb$query)
  }
  gb
}
