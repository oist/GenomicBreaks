#' Tamura 1992 distance
#'
#' The Tamura 1992 (T92) distance extends the K80 distance by taking `GC` content into account.  It is calculated as \eqn{-h \ln \left(1 - \frac{p}{h} - q\right) - \frac{1}{2} \times (1 - h) \ln\left(1 - 2 q\right)}, where \eqn{p} is the probability of transition, \eqn{q} the probability of transversion, \eqn{h = 2\theta (1 - \theta)} and \eqn{\theta} is the `GC` content. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#T92_model_(Tamura_1992)) for more details.
#'
#' @param train_parameters A list containing the probabilities of the alignment, produced by the [`readTrainFile()`] function.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#'
#' @family Similarity indexes
#'
#' @author Zikun Yang
#'
#' @references Tamura, K. (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases." *Molecular Biology and Evolution*, 9(4), 678–687. DOI: [10.1093/oxfordjournals.molbev.a040752](https://doi.org/10.1093/oxfordjournals.molbev.a040752)#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' T92_distance(parameters)
#'
#' @export

T92_distance <- function(train_parameters) {
  theta <- train_parameters[["probability_G_G"]] + train_parameters[["probability_C_C"]]
  h <- 2 * theta * (1 - theta)
  p <- train_parameters[["probability_A_G"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_C_T"]]
  q <- train_parameters[["probability_A_C"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_G_C"]]
  - h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q)
}
