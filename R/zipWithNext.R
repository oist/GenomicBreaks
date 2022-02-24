#' Zip elements with the next ones
#'
#' Zip an object with itself after inserting a variable number of steps.  For
#' instance, with a step of 1, _a_ _b_ _c_ would become _a,b_, _b,c_, _c,NA_.
#'
#' This is a programming pattern that I use a lot in this package, to prepare
#' vectors fro parallel processing.  I typically use it on metadata columns
#' of `GBreaks` objects.  I created this function because I always forget if I
#' should use `head` or `tail`, which causes bugs or delays in development.
#'
#' @param x Object with a vector structure.
#' @param step Remove the first `step` entries from the object before zipping
#'        it with itself.
#'
#' @return A [`S4Vectors::Pairs`] object with a copy of `x` in the `first` slot,
#' and `x` with `step` elements removed from its head and complemented with
#' `NAs` in the `second` slot.
#'
#'
#' @examples
#' zipWithNext(LETTERS, step = 2)
#'
#' @importFrom S4Vectors Pairs
#' @importFrom utils head tail
#' @export

zipWithNext <- function(x, step = 1) {
  Pairs(
    x,
    c(tail(x, step * -1), rep(NA, step))
  )
}
