#' The GenomicBreaks class
#'
#' The `GBreaks` class is a simple wrapper to the `GRanges` class
#'
#' @export

GBreaks <- setClass("GBreaks", contains = "GRanges")
