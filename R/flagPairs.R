# Table of all combinations
allPossiblePairClasses <- expand.grid(
  Qnext = c("next", "prev", "other"),
  Strnd = c("+", "-"),
  Snext = c("+", "-")
)
allPossiblePairClasses$pair <- NA
allPossiblePairClasses <- with(allPossiblePairClasses, {
  # Colinearity
  allPossiblePairClasses$pair[Qnext == "next" & Strnd == "+" & Snext == "+"] <- "Col"
  allPossiblePairClasses$pair[Qnext == "prev" & Strnd == "-" & Snext == "-"] <- "Col"
  # Inversion
  allPossiblePairClasses$pair[Qnext == "next" & Strnd == "+" & Snext == "-"] <- "Inv"
  allPossiblePairClasses$pair[Qnext == "next" & Strnd == "-" & Snext == "+"] <- "Inv"
  allPossiblePairClasses$pair[Qnext == "prev" & Strnd == "+" & Snext == "-"] <- "Inv"
  allPossiblePairClasses$pair[Qnext == "prev" & Strnd == "-" & Snext == "+"] <- "Inv"
  # "Flipped", like A jumped in front of B
  allPossiblePairClasses$pair[Qnext == "prev" & Strnd == "+" & Snext == "+"] <- "Flp"
  allPossiblePairClasses$pair[Qnext == "next" & Strnd == "-" & Snext == "-"] <- "Flp"
  # "Scrambled" when query ranges are not next to each other.
  allPossiblePairClasses$pair[Qnext == "other"] <- "Scr"
  allPossiblePairClasses$paste <- paste(Qnext, Strnd, Snext)
  allPossiblePairClasses
})

#' Flag successive entries of a `GBreaks` object
#'
#' Scans a sorted [`GBreaks`] object and issues a flag describing the relation
#' between the current entry and the next entry.
#'
#' Because the object is sorted, the next entry on the _target_ genome is by
#' definition following the current one unless the end of the sequence feature
#' (contig, scaffold, …) is reached.  The two _target_ ranges can be on the
#' same strand or on opposite strands.  No difference is made between `+/-`
#' and `-/+` orientations because it will be arbitrary unless the sequence
#' features of both genomes are oriented, as in the case of comparing
#' chromosomal assemblies of related species.
#'
#' Depending on whether the _query_ ranges precede or follow each other, or
#' are not next to each other, and depending on the strand of the alignments,
#' a flag is issued on the current entry.  `Col` signals colinearity with the
#' next entry, `Inv` an inversion on either entry, `Flp` signals that the order
#' of the entries is as if one hand jumped over the other one.  When the _query_
#' ranges are not next to each other, a `Scr` flag, for _scrambled_, is issued.
#' Lastly, `Bnd` (for _boundary_) signals that there is no pair to analyse
#' because the current entry is the last one for the current sequence level on
#' the _target_ genome.
#'
#' @param gb A `GBreaks()` object.
#'
#' @return Returns the `GBreaks` object with an extra `pairs` metadata column.
#' This `pairs` column is a factor of all flags, so that `table(gb$pairs)`
#' always returns a vector of the same length, reporting `0` when necessary.
#' This allows to easily aggregate results from multiple objects.
#'
#' @family Flagging functions
#' @family Inversion functions
#' @family Colinearity functions
#' @family Structural variants
#'
#' @examples
#' flagPairs(exampleInversion)
#' flagPairs(exampleTranslocation)
#'
#' # How the pair classes are defined:
#' GenomicBreaks:::allPossiblePairClasses
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede follow
#' @export

flagPairs <- function (gb) {
  if(length(gb) == 0) {
    gb$pairs <- factor(character(0))
    return(gb)
  }
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")
  gb.bak <- gb # save the original object

  # Compute various indicators
  gb$Strnd <- as.character(strand(gb))
  gb$Snext <- c(tail(as.character(strand(gb)), -1), NA)
  gb$oneNext <- precede(gb$query, tail(gb$query, -1), ignore.strand = T) - seq_along(gb) == 0
  gb$onePrev <- follow( gb$query, tail(gb$query, -1), ignore.strand = T) - seq_along(gb) == 0
  gb$Qnext <- "other"
  gb$Qnext[sapply(gb$oneNext, isTRUE)] <- "next"
  gb$Qnext[sapply(gb$onePrev, isTRUE)] <- "prev"

  # Match the combinations with the allPossiblePairClasses table
  gb$paste <- paste(gb$Qnext, gb$Strnd, gb$Snext)
  gb$paste[is.na(precede(gb, ignore.strand=T))] <- NA   # Put NAs at boundaries
  gb$paste <- factor(gb$paste, labels = allPossiblePairClasses$pair, levels = allPossiblePairClasses$paste)

  # Transform the NAs into the "boundary" state
  gb$paste <- as.character(gb$paste)
  gb$paste[is.na(gb$paste)] <- "Bnd"

  # Output as a factor
  gb.bak$pairs <- factor(gb$paste, levels = c("Col", "Inv", "Flp", "Scr", "Bnd"))
  gb.bak
}
