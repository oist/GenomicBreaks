#' Tamura 1992 distance
#'
#' The Tamura 1992 (T92) distance extends the K80 distance by taking `GC` content into account.  It is calculated as \eqn{-h \ln \left(1 - \frac{p}{h} - q\right) - \frac{1}{2} \times (1 - h) \ln\left(1 - 2 q\right)}, where \eqn{p} is the probability of transition, \eqn{q} the probability of transversion, \eqn{h = 2\theta (1 - \theta)} and \eqn{\theta} is the `GC` content. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#T92_model_(Tamura_1992)) for more details.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @family Nucleotide distances
#'
#' @author Zikun Yang
#'
#' @references Tamura, K. (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases." *Molecular Biology and Evolution*, 9(4), 678–687. DOI: [10.1093/oxfordjournals.molbev.a040752](https://doi.org/10.1093/oxfordjournals.molbev.a040752)
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' T92_distance(exampleSubstitutionMatrix)
#'
#' @export

T92_distance <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  m <- m[c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T')]
  P <- prop.table(m)
  theta <- P['G', 'G'] + P['C', 'C']
  h <- 2 * theta * (1 - theta)
  p <- P['A', 'G'] + P['G', 'A'] +
       P['C', 'T'] + P['T', 'C']
  q <- P['A', 'C'] + P['C', 'A'] +
       P['G', 'T'] + P['T', 'G'] +
       P['A', 'T'] + P['T', 'A'] +
       P['C', 'G'] + P['G', 'C']
  - h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q)
}
