#' Equilibrium GC fraction (GC*)
#'
#' Computes the predicted equilibrium GC content \strong{GC*} from a nucleotide
#' substitution count matrix.  GC* is the fraction of G/C toward which base
#' composition will evolve under a two-class (weak vs. strong) mutation-bias
#' model, using the asymmetry between AT→GC (W→S) and GC→AT (S→W) changes
#' (Sueoka, 1962).
#'
#' The classic directional-mutation equilibrium for GC content is:
#'
#' \deqn{\mathrm{GC}^* = \frac{\text{W}\to\text{S}}{\text{W}\to\text{S} \quad+\quad \text{S}\to\text{W}}}{GC* = W→S / (W→S + S→W)}
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @note If the _target_ genome is not the true ancestor (which is likely in
#' simple pairwise comparisons of extant genomes), GC* should be interpreted
#' cautiously.  It does not predict future GC content, but its position relative
#' to the GC of the _target_ can still indicate the direction of substitution
#' bias: if GC* is higher than the target GC, the bias favors G/C; if lower,
#' it favors A/T.  Transpose the input matrix to study the _query_ genome the
#' same way.
#'
#' @references Noboru Sueoka.  On the genetic basis of variation and
#' heterogeneity of DNA base composition. (1962) *Proc Natl Acad Sci U S A*
#' 48(4):582-92. \doi{10.1073/pnas.48.4.582}
#'
#' @returns A numeric value between 0 and 1.
#'
#' @family Nucleotide distances
#'
#' @author Charles Plessy
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' GCequilibrium(parameters)
#' GCequilibrium(t(parameters$probability_matrix))
#'
#' @export

GCequilibrium <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  W <- c("A", "T")
  S <- c("G", "C")
  a <- sum(m[W, S])
  b <- sum(m[S, W])
  a / (a + b)
}
