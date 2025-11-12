#' Hasegawa–Kishino–Yano (1985) distance (HKY85)
#'
#' Computes the HKY85 nucleotide substitution distance between two sequences,
#' using only the A/C/G/T block of a pairwise count matrix. Base frequencies
#' are estimated as the average of row and column marginals, and transitions
#' are pooled into a single class (A↔G and C↔T).
#'
#' @details
#' Let \eqn{P} be the total fraction of **transition** differences
#' (A↔G **plus** C↔T) among A/C/G/T pairs, and \eqn{Q} the fraction of
#' **transversion** differences (all other A↔\{C,T\}, G↔\{C,T\} pairs).
#' Let \eqn{\pi_A,\pi_C,\pi_G,\pi_T} be the **average base frequencies**
#' across the two sequences (row and column marginals averaged),
#' \eqn{\pi_R = \pi_A + \pi_G}, and \eqn{\pi_Y = \pi_C + \pi_T}.
#'
#' The HKY85 distance (special case of TN93 with a single transition rate) is:
#' \deqn{
#' d_{\mathrm{HKY85}} =
#'   -\,2(\pi_A \pi_G + \pi_C \pi_T)\;
#'      \ln\!\left(1 - \frac{P}{2(\pi_A \pi_G + \pi_C \pi_T)}\right)
#'   -\,2\,\pi_R \pi_Y\;
#'      \ln\!\left(1 - \frac{Q}{2\,\pi_R \pi_Y}\right).
#' }
#'
#' Implementation notes:
#' \itemize{
#'   \item Rows/columns outside \code{A,C,G,T} are dropped.
#'   \item \eqn{P} and \eqn{Q} are computed from the normalized 4×4 table.
#'   \item \eqn{\pi_\bullet} are estimated as the average of row and column marginals.
#'   \item Returns \code{Inf} if any log argument is non‑positive or a denominator is zero.
#' }
#'
#' @param m A numeric matrix of pair counts with row/column names including
#'   \code{A,C,G,T}. If \code{m} is a list, the field \code{$probability_matrix}
#'   is used (for compatibility with your other helpers).
#'
#' @return A single numeric: the HKY85 distance (substitutions/site).
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @references
#' Hasegawa, M., Kishino, H., and Yano, T. (1985). Dating of the human–ape
#' splitting by a molecular clock of mitochondrial DNA. *J. Mol. Evol.* 22:160–174.
#'
#' RevBayes tutorial (overview of HKY, parameters \eqn{\kappa,\pi}):
#' \url{https://revbayes.github.io/tutorials/ctmc/}
#'
#' IQ‑TREE model manual (relationships among HKY, K80, TN93):
#' \url{https://iqtree.github.io/doc/Substitution-Models}
#'
#' @author Entirely written by Copilot and not proofchecked yet
#'
#' @examples
#' HKY85_distance(exampleSubstitutionMatrix)
#'
#' @export

HKY85_distance <- function(m) {
  if (is.list(m)) m <- m$probability_matrix
  m <- m[c("A", "C", "G", "T"), c("A", "C", "G", "T")]
  if (all(m == 0)) return(NA)
  P <- prop.table(m)

  # Average base frequencies across the pair
  target_freqs <- rowSums(P)
  query_freqs  <- colSums(P)
  gA <- (target_freqs["A"] + query_freqs["A"]) / 2
  gC <- (target_freqs["C"] + query_freqs["C"]) / 2
  gG <- (target_freqs["G"] + query_freqs["G"]) / 2
  gT <- (target_freqs["T"] + query_freqs["T"]) / 2
  gR <- gA + gG
  gY <- gC + gT

  # Observed proportions
  PR <- P["A", "G"] + P["G", "A"]           # purine transitions
  PY <- P["C", "T"] + P["T", "C"]           # pyrimidine transitions
  P_tr <- PR + PY                           # total transitions (HKY has one transition rate)
  Q    <- P["A", "C"] + P["C", "A"] +       # transversions
    P["A", "T"] + P["T", "A"] +
    P["G", "C"] + P["C", "G"] +
    P["G", "T"] + P["T", "G"]

  # HKY denominators (specialization of TN93 with PR & PY combined)
  S <- gA * gG + gC * gT        # π_A π_G + π_C π_T
  a <- 2 * S
  c <- 2 * gR * gY

  t1 <- 1 - P_tr / a
  t2 <- 1 - Q     / c

  # Guard against invalid logs / zero denominators
  if (a <= 0 || c <= 0 || t1 <= 0 || t2 <= 0) return(Inf)

  unname(-a * log(t1) - c * log(t2))
}
