#' Guesstimate seqlevel lengths
#'
#' When no [`seqlengths`] are available, one can resort to set them as to the
#' maximal end coordinate found in the object.
#'
#' @param gr A [`GRanges`] object
#'
#' @returns [`guessSeqLengths`] returns sequence lengths that have been guessed
#' as the maximal end coordinates found in the object, or the object's
#' `seqlengths` if if they did already exist.
#'
#' @author Charles Plessy
#'
#' @family modifier functions
#'
#' @examples
#' exampleTranslocation
#' seqlengths(exampleTranslocation$query)
#' guessSeqLengths(exampleTranslocation$query)
#' gb2 <- exampleTranslocation
#' seqlengths(gb2$query) <- NA
#' guessSeqLengths(gb2$query)
#'
#' @export

guessSeqLengths <- function(gr) {
  if (!any(is.na(seqlengths(gr)))) return (seqlengths(gr))
  tapply(end(gr), seqnames(gr), max)
}

#' @rdname guessSeqLengths
#'
#' @returns [`forceSeqLengths`] returns the object in which the sequence lengths
#' have been set to the maximal end coordinates found in the object if if they
#' did not exist.
#'
#' @family modifier functions
#'
#' @examples
#' forceSeqLengths(exampleTranslocation$query)
#' forceSeqLengths(exampleTranslocation$query) |> seqlengths()
#'
#' @export

forceSeqLengths <- function(gr) {
  if (!any(is.na(seqlengths(gr)))) return (gr)
  seqlengths(gr) <- as.vector(guessSeqLengths(gr))
  gr
}
