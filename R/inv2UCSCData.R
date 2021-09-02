#' Represent inversions as UCSCData objects
#'
#' Scan a [`GBreaks`] object for inversions and return UCSC track object
#' in [`rtracklayer::UCSCData-class`] format for export in BED12 format with
#' a command such as `rtracklayer::export(x, "test.bed", "BED")`.
#'
#' @param gb A `GBreaks` object.
#'
#' @returns Returns a `UCSCData` object.
#'
#' @export

inv2UCSCData <- function(gb) {
  gb <- flagInversions(gb)
  if (sum(gb$inv) == 0) stop("No inversion found!")

  # Index of the first range in each inversion triplet.
  invPos <- which(gb$inv)

  # Function to build one representation for one inversion
  triplet2line <- function(index, gb) {
    triplet <- c(index, index + 1, index + 2)
    gb[triplet]
    invObj <- range(gb[triplet], ignore.strand = TRUE)
    strand(invObj) <- strand(gb[triplet][2])
    invObj <- as(invObj, "UCSCData")
    invObj$blocks <-
      IRangesList(ranges(gb[triplet]) |> shift(1 - min(start(gb[triplet]))))
    invObj$thick <- ranges(gb[triplet][2])
    invObj@trackLine <- new( "BasicTrackLine", name = "Inversions"
                           , description = "GenomicBreaks inversions"
                           , visibility="full")
    names(invObj) <- as.character(gb[triplet][2])
    invObj$itemRgb <- "black"
    invObj
  }

  # Apply the function to all inversions and aggregate results.
  do.call(c, lapply(invPos, triplet2line, gb))
}
