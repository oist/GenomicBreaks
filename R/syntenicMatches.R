#' Sequence levels matched by synteny
#'
#' @param gb A [`GBreaks`] object.
#' @param min A proportion threshold between 0 and 1 for a pair of sequence
#'        levels to be considered matching.
#'
#' @family Synteny functions
#'
#' @returns `syntenicProportions` returns a list named by the sequence levels
#' of the _target_ genome, indicating what proportion of their aligned sequences
#' are mapped to a given sequence level of the _query_ genome.
#' `syntenicMatches` only retains the names of the _query_ sequence levels after
#' discarding the ones representing a proportion of hits lower than the `min`
#' threshold.
#'
#' @examples
#' syntenicProportions(exampleTranslocation)
#' syntenicMatches(exampleTranslocation)
#'
#' @export

#' @rdname syntenicMatches

syntenicProportions <- function(gb) {
  if(length(gb) == 0) return(list())
  gbl <- split(gb, seqnames(gb), drop = TRUE)
  l <- lapply(names(gbl), \(x) {
    x <- gbl[[x]]
    proportions <- (split(x, seqnames(x$query)) |> width() |> sum() |> sort(decreasing = TRUE)) / sum(width(x))
  })
  names(l) <- names(gbl)
  l
}

#' @rdname syntenicMatches

syntenicMatches <- function(gb, min=0.01) {
  l <- syntenicProportions(gb)
  lapply(l, \(x) {
    names(x[x >= min])
  })
}
