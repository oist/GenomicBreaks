#' Load one-to-one orthogroups for a set of species
#'
#' Takes a _phylogenetic hierarchical orthogroup_ table from the output of
#' _OrthoFinder_, optionally subsets it for some species, and returns the
#' table containing one-to-one orthogroups.
#'
#' Orthogroups that are not one-to-one have either a missing entry for one
#' species, represented as an empty string, or multiple entries, separated
#' with comma characters.  Species that are not part of the clade covered by
#' the table have `NA` in every row.
#'
#' @param file Path to the _OrthoFinder_ table.
#'
#' @param species A character vector of species annotation names matching
#' exactly the names of the columns in the _OrthoFinder_ table.
#'
#' @returns Returns a [`DataFrame`] with one column per species, containing
#' protein identifiers, and the columns `HOG`, `OG` and `Gene.Tree.Parent.Clade`.
#'
#' @export

load_one_to_ones <- function(file, species = NULL) {
  # Read table
  pho <- read.delim(file, check.names = FALSE)

  # Remove all-NA columns
  pho <- pho[,sapply(pho, \(col) ! all(is.na(col)))]

  # Subset for species list
  if (! is.null(species))
    pho <- pho[,c("HOG", "OG", "Gene Tree Parent Clade", species)]

  # Find unique entries
  # An entry is unique if no comma character is found in
  un <- do.call(cbind, lapply(pho,\(x) ! grepl(",", x)))

  # Find non-empty entries
  nonempt <- do.call(cbind, lapply(pho, \(x) x != ""))

  # Flag rows containing only unique entries
  un.rows <- apply(un & nonempt, 1, all)

  # Remove non-unique entries
  pho.un <- pho[un.rows,]

  # Return as a DataFrame
  as(pho.un, "DataFrame")
}
