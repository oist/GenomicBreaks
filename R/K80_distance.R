#' Index representing the evolutionary distance between two genomes.
#'
#' The K80 distance is calculated as -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q)). See [WIKI](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)) for more details.
#'
#' @param train_parameters A list containing the probabilities of the alignment.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.ss
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' K80_distance(parameters)
#'
#' @export

# K80 model (Kimura 1980)
K80_distance <- function(train_parameters) {

  p <- train_parameters[["probability_A_G"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_C_T"]]
  q <- train_parameters[["probability_A_C"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_G_C"]]
  k <- - 0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
  return(k)

}
