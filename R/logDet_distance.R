#' LogDet (log-determinant) distance
#'
#' @description
#' Computes the LogDet (a.k.a. log-determinant) nucleotide distance between two
#' sequences from a 4×4 A/C/G/T joint count table. The method normalizes the
#' joint frequency matrix by the marginal base compositions of the two sequences
#' and then applies the log-determinant transform.
#'
#' @details
#' Let \eqn{M} be the 4×4 **joint count** matrix over \code{A,C,G,T} and
#' \eqn{F = M / \sum M} the **joint frequency** (divergence) matrix. Let
#' \eqn{\boldsymbol{r} = F \mathbf{1}} (row sums, seq1 base frequencies) and
#' \eqn{\boldsymbol{c} = F^\top \mathbf{1}} (column sums, seq2 base frequencies).
#' Define diagonal matrices \eqn{D_R = \mathrm{diag}(\boldsymbol{r})} and
#' \eqn{D_C = \mathrm{diag}(\boldsymbol{c})}. The **LogDet** normalization uses
#' \deqn{
#' S \;=\; D_R^{-1/2}\, F \, D_C^{-1/2},
#' }
#' and the **distance** is
#' \deqn{
#' d_{\mathrm{LogDet}} \;=\; -\,\ln \!\big( \det S \big)
#' \;=\; -\,\ln \!\big( \det F \big) \;+\; \tfrac{1}{2}\,\big[ \ln \!\big( \det D_R \big)
#' \;+\; \ln \!\big( \det D_C \big) \big].
#' }
#' This formulation follows the LogDet/paralinear family where determinants
#' commute (real scalars), enabling robustness to composition heterogeneity; see
#' Lockhart et al. (1994) for the transform and discussion, and MEGA’s note on
#' non-computability when log terms approach zero. [1](https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf)[2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)
#'
#' @param m A numeric matrix of pairwise counts with row/column names including
#'   \code{A,C,G,T}. If \code{m} is a list, the field \code{$probability_matrix}
#'   is used (compatible with your other helpers).
#' @param pseudocount A non-negative number added to each A/C/G/T cell before
#'   normalization (default \code{0}). Use a small value (e.g. \code{0.5}) to
#'   reduce the chance of a singular matrix when some joint cells are zero. [2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)
#'
#' @returns A single numeric: the LogDet distance (substitutions/site, unitless).
#'   Returns \code{Inf} if any required determinant/sign is non‑positive.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @return A single numeric: the HKY85 distance (substitutions/site).
#'
#' @references
#' Lockhart PJ, Steel MA, Hendy MD, Penny D (1994). Recovering evolutionary
#' trees under a more realistic model of sequence evolution. *Mol Biol Evol*
#' **11**:605–612. \url{https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf} [1](https://evoluscope.fr/phylographe/biblio/Lockhartetal1994.pdf)
#'
#' MEGA help: “LogDet Distance Could Not Be Computed” (notes on log terms).
#' \url{https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm} [2](https://www.megasoftware.net/web_help_11/logDet_distance_Could_Not_Be_Computed.htm)
#'
#' Gu X, Li WH (1996). Bias-corrected paralinear and LogDet distances and
#' tests under nonstationary nucleotide frequencies. *Mol Biol Evol* **13**:1375–1383.
#' \url{https://paperity.org/p/54207439/bias-corrected-paralinear-and-logdet-distances-and-tests-of-molecular-clocks-and} [3](https://paperity.org/p/54207439/bias-corrected-paralinear-and-logdet-distances-and-tests-of-molecular-clocks-and)
#'
#' @examples
#' logDet_distance(exampleSubstitutionMatrix)
#'
#' @export

logDet_distance <- function(m, pseudocount = 0) {
  if (is.list(m)) m <- m$probability_matrix
  # Keep only A/C/G/T
  nts <- c("A", "C", "G", "T")
  m <- m[nts, nts, drop = FALSE]

  # Optional smoothing to avoid zero cells
  if (pseudocount > 0) {
    m <- m + pseudocount
  }

  # Joint frequency (divergence) matrix
  total <- sum(m)
  if (total <= 0) return(Inf)
  F <- m / total

  # Marginal base compositions
  r <- rowSums(F)
  c <- colSums(F)

  # Guard: all marginals must be positive
  if (any(r <= 0) || any(c <= 0)) return(Inf)

  # Normalization: S = D_R^{-1/2} F D_C^{-1/2}
  inv_sqrt_r <- 1 / sqrt(r)
  inv_sqrt_c <- 1 / sqrt(c)
  S <- (inv_sqrt_r * F) * rep(inv_sqrt_c, each = 4)  # fast diag scaling

  # Compute log(det(S)) robustly
  detS <- determinant(S, logarithm = TRUE)
  if (detS$sign <= 0) return(Inf)
  logdetS <- as.numeric(detS$modulus)

  # Distance
  -logdetS
}
