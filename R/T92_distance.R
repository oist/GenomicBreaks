#' Index representing the evolutionary distance between two genomes.
#'
#' The T92 distance is calculated as -h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q). See [WIKI](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#T92_model_(Tamura_1992)) for more details.
#' @param train_parameters A list containing the probabilities of the alignment.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' T92_distance(parameters)
#'
#' @export

# T92 model (Tamura 1992)
T92_distance <- function(train_parameters) {

  theta <- train_parameters[["probability_G_G"]] + train_parameters[["probability_C_C"]]
  h <- 2 * theta * (1 - theta)
  p <- train_parameters[["probability_A_G"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_C_T"]]
  q <- train_parameters[["probability_A_C"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_G_C"]]
  d <- - h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q)
  return(d)

}
