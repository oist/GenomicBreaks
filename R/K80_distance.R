#' Kimura's 2-parameter distance
#'
#' The Kimura's 2-parameter (K80) distance distinguishes between transitions and transversions.  It is calculated as \eqn{-\frac{1}{2} \ln\left(\left(1 - 2p - q\right) \times \sqrt{1 - 2q}\right)}, where \eqn{p} is the probability of transition and \eqn{q} is the probability of transversion. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)) for more details.
#'
#' @param train_parameters A list containing the probabilities of the alignment, produced by the [`readTrainFile()`] function.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#'
#' @family Similarity indexes
#'
#' @author Zikun Yang
#'
#' @references Kimura, M. (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." *Journal of Molecular Evolution*, 16, 111–120. DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581)
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' K80_distance(parameters)
#'
#' @export

K80_distance <- function(train_parameters) {
  p <- train_parameters[["probability_A_G"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_C_T"]]
  q <- train_parameters[["probability_A_C"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_G_C"]]
  - 0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
}
