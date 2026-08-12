#' Inversion Distance
#'
#' Computes the minimal number of inversions required to sort a signed permutation
#' using the Hannenhalli and Pevzner algorithm.
#'
#' This function uses several helper functions defined in breakpointGraph.R (e.g., `extendedPermutation`,
#' `breakpoint_graph`, `hurdles_count`, `superhurdles_count`, and others) to
#' compute properties of the breakpoint graph and identify cycles, hurdles and superhurdles.
#' It also depends on `permutationVector()`, which is defined in permutationVector.R.
#'
#' This algorithm was designed to work in *a single, linear chromosome alignment*.
#' Although the function still works if the GBreaks object involves more than one chromosome,
#' the returned value for the minimal number of inversions will imply in non-usual inversions
#' if different chromosomes have orthologous regions.
#'
#' @references Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into turnip: polynomial algorithm for sorting signed permutations by reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.
#'
#' @param x Either a GBreaks object or a signed permutation vector.
#'   If \code{x} is a GBreaks object, a permutation vector will be extracted
#'   using \code{\link{permutationVector}}.
#'
#' @return An integer: the minimal number of inversions needed to sort the permutation.
#'
#' @examples
#' \dontrun{
#' # Example using a permutation vector directly
#' # Suppose we want to sort the permutation p = c(1, 3, -2, 4)
#' inversionDistance(c(1, 3, -2, 4))
#'
#' # Example using a GBreaks object
#' example         <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690"))
#' strand(example) <- c(              "+",            "-",            "-",            "-",            "+",            "+"      )
#' example$query   <- GRanges(c("chrA:100-190", "chrA:300-390", "chrA:200-290", "chrA:400-490", "chrA:600-690", "chrA:500-590"))
#' gb_example      <- GBreaks(example)
#' inversionDistance(gb_example)}
#'
#' @seealso \code{\link{permutationVector}} for generating the permutation vector.
#'
#' @family Rearrangement distances
#' @family Similarity indexes
#'
#' @importFrom igraph make_empty_graph add_edges E V E<- V<- components induced_subgraph ends delete_vertices
#'
#' @author Bruna Fistarol
#'
#' @export
inversionDistance <- function(x){

  if (inherits(x, "GBreaks")) {
    p <- permutationVector(x)

  } else if (is.numeric(x) && is.vector(x)) {
    if (any(x != as.integer(x))) {
      stop("Permutation vector must contain only integers.")
    }
    p <- as.integer(x)

    n <- length(p)
    if (!identical(sort(abs(p)), 1:n)) {
      stop("Permutation vector must contain each of 1 to n exactly once (up to sign).")
    }

  } else {
    stop("Input must be either a GBreaks object or a signed permutation vector.")
  }

  if (identical(p, seq_len(length(p)))) {
    return(0L)
  }

  p_extended <- extendedPermutation(p)

  bp_graph <- breakpoint_graph(p_extended)
  hurdles <- hurdles_count(bp_graph, p_extended)
  superhurdles <- superhurdles_count(hurdles, bp_graph, p_extended)

  bp_count(p_extended) - cycle_count(bp_graph) + sum(hurdles$hurdle) + is_fortress(superhurdles)

}
