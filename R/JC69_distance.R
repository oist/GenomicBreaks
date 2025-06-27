#' Index representing the evolutionary distance between two genomes.
#'
#' The JC69 distance is calculated as -0.75 * log(1 - 4 * p / 3), where p is the sum of the probabilities of the alignment. See [WIKI](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) for more details.
#'
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
#' JC69_distance(parameters)
#'
#' @export

# JC69 model (Jukes and Cantor 1969)
JC69_distance <- function(train_parameters) {

  p <- train_parameters[["probability_A_C"]] + train_parameters[["probability_A_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_C_T"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_G_C"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_T_G"]]
  d <- - 0.75 * log(1 - 4 * p / 3)
  return(d)

}
