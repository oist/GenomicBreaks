#' Bridge regions
#'
#' Maps unaligned regions of the _target_ to the _query_ genome when they are
#' flanked by colinear regions.
#'
#' @note Because some aligned regions can be directly adjacent (no gaps), the
#' returned `GBreaks` object may contain ranges of width zero, where the _start_
#' coordinate is 1 nucleotide higher than the _end_ coordinate.
#'
#' @references Bridge regions have also been called \dQuote{simultaneous gaps}
#' in the comparison of the mouse and human genomes by Kent WJ, Baertsch R,
#' Hinrichs A, Miller W, Haussler D. (_Evolution's cauldron: duplication,
#' deletion, and rearrangement in the mouse and human genomes._ Proc Natl Acad
#' Sci U S A. 2003;100(20):11484-11489. doi:10.1073/pnas.1932072100)
#'
#' @param gb A [`GBreaks`] object.
#'
#' @return Returns a new `GBreaks` object of shorter length.
#'
#' @family Colinearity functions
#' @family modifier functions
#'
#' @author Charles Plessy
#'
#' @examples
#' exampleColinear5
#' bridgeRegions(exampleColinear5)
#'
#' # Note the zero-width ranges when aligned regions are directly adjacent.
#' exampleColinear3
#' bridgeRegions(exampleColinear3)
#'
#' @export

bridgeRegions <- function(gb) {
  # Collect colinear regions and discard the rest
  colinearRegions <- filterColinearRegions(flagColinearAlignments(gb), rename = FALSE)
  # Turn the runs of [(TRUE)n, FALSE]n into indices identifying each colinear region
  idx <- c(0, head(cumsum(!colinearRegions$colinear), -1))
  # Split into a GRangesList
  gbl <- split(colinearRegions, idx)
  # Collect gap ranges in the target and query genomes of each colinear regions.
  br <- endoapply(gbl, \(gb) {
    # Check strand
    onMinus <- all(strand(gb) == "-")
    # Need to subtract 1 temporarly because some ranges are adjacent (no gap).
    GBreaks(target =                            cleanGaps(gb      - 1) -1 ,
            query  = sort(decreasing = onMinus, cleanGaps(gb$query -1) -1))
  }) |> unlist()
  # Remove names and return
  names(br) <- NULL
  br
}
