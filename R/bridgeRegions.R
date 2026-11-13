#' Bridge regions
#'
#' Maps unaligned regions of the _target_ to the _query_ genome when they are
#' flanked by colinear regions.
#'
#' @note Because some aligned regions can be directly adjacent (no gaps), the
#' returned `GBreaks` object may contain ranges of width zero, where the _start_
#' coordinate is 1 nucleotide higher than the _end_ coordinate.
#'
#' @references The term \dQuote{bridge region} was introduced in: Plessy C,
#' Mansfield MJ, Bliznina A, Masunaga A, West C, Tan Y, Liu AW, Grašič J, Del
#' Río Pisula MS, Sánchez-Serna G, Fabrega-Torrus M, Ferrández-Roldán A,
#' Roncalli V, Navratilova P, Thompson EM, Onuma T, Nishida H, Cañestro C,
#' Luscombe NM (2024). Extreme genome scrambling in marine planktonic
#' _Oikopleura dioica_ cryptic species. *Genome Research*, 34(3), 426–440.
#' \url{https://doi.org/10.1101/gr.278295.123}. PMID: 38621828
#'
#' Bridge regions have also been called \dQuote{simultaneous gaps} in the
#' comparison of the mouse and human genomes by Kent WJ, Baertsch R, Hinrichs
#' A, Miller W, Haussler D. (_Evolution's cauldron: duplication, deletion, and
#' rearrangement in the mouse and human genomes._ Proc Natl Acad Sci U S A.
#' 2003;100(20):11484-11489. doi:10.1073/pnas.1932072100)
#'
#' @param gb A [`GBreaks`] object.
#'
#' @return Returns a new `GBreaks` object of shorter length.  Its length will
#' be zero if no bridge regions were found.
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
  if (length(colinearRegions) == 0) return (colinearRegions)
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
