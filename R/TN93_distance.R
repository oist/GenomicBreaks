#' Tamura–Nei 1993 distance
#'
#' The Tamura–Nei 1993 (TN93) distance extends the T92 distance by allowing
#' unequal base frequencies and two distinct transition components, plus
#' transversions.
#'
#' Let \eqn{P_R} be the fraction of purine transitions (A↔G), \eqn{P_Y} the
#' fraction of pyrimidine transitions (C↔T), and \eqn{Q} the fraction of
#' transversions among A/C/G/T pairs.  Let \eqn{\pi_A, \pi_C, \pi_G, \pi_T} be
#' the average base frequencies across the two sequences (row and column
#' marginals averaged), \eqn{\pi_R = \pi_A + \pi_G}, and \eqn{\pi_Y = \pi_C + \pi_T}.
#'
#' The TN93 distance is:
#' \deqn{
#' d_{\mathrm{TN93}} =
#'   -\,2 \pi_A \pi_G \, \ln\!\left(1 - \frac{P_R}{2 \pi_A \pi_G}\right)
#'   -\,2 \pi_C \pi_T \, \ln\!\left(1 - \frac{P_Y}{2 \pi_C \pi_T}\right)
#'   -\,2 \pi_R \pi_Y \, \ln\!\left(1 - \frac{Q}{2 \pi_R \pi_Y}\right)
#' }
#'v
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' aligned to bases on the _query_ genome.  As a convenience it can also receive
#' a list produced by the [`readTrainFile()`] function, containing this matrix.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @author Charles Plessy / M365 Copilot (GPT-5 on)
#'
#' @references
#' Tamura, K. and Nei, M. (1993).  Estimation of the number of nucleotide
#' substitutions in the control region of mitochondrial DNA in humans and
#' chimpanzees.  *Molecular Biology and Evolution* **10**:512–26.
#' \doi{10.1093/oxfordjournals.molbev.a040023}
#'
#' @returns Returns a numeric value show the evolutionary distance between two
#' genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' TN93_distance(exampleSubstitutionMatrix)
#'
#' @export

TN93_distance <- function(m) {
  if (is.list(m)) m <- m$probability_matrix
  m <- m[c("A", "C", "G", "T"), c("A", "C", "G", "T")]
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
  PR <- P["A", "G"] + P["G", "A"]           # transitions within purines
  PY <- P["C", "T"] + P["T", "C"]           # transitions within pyrimidines
  Q  <- P["A", "C"] + P["C", "A"] +         # transversions
        P["A", "T"] + P["T", "A"] +
        P["G", "C"] + P["C", "G"] +
        P["G", "T"] + P["T", "G"]

  # TN93 denominators
  a <- 2 * gA * gG
  b <- 2 * gC * gT
  c <- 2 * gR * gY

  t1 <- 1 - PR / a
  t2 <- 1 - PY / b
  t3 <- 1 - Q  / c

  # Guard against invalid logs / zero denominators
  if (a <= 0 || b <= 0 || c <= 0 || t1 <= 0 || t2 <= 0 || t3 <= 0) return(Inf)

  unname(-a * log(t1) - b * log(t2) - c * log(t3))
}
