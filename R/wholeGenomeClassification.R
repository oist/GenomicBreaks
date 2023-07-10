#' Whole-genome classification object
#'
#' Classifies aligned genomic regions as _isolated_ or _collinear_, and
#' unaligned regions as _breakpoint_ or _bridge_ regions.  The _bridge_ regions
#' are always flanked by _collinear alignments_ and the _isolated alignments_
#' are always flanked by _breakpoint regions_.
#'
#' @param gb A [`GBreaks`] object representing a one-to-one whole genome
#'        alignment.
#' @param coa The coalesced one-to-one alignment.  If not provided, it will be\
#'        computed on-the-fly with the [`coalesce_contigs`] function.
#'
#' @returns A [`GenomicRanges::GRanges`] object representing the _target_
#' genome, inheriting its _sequence information_ ([`GenomeInfoDb::Seqinfo`]).
#' The class of each region is indicated by a factor in the `type` metadata
#' column.
#'
#' @author Charles Plessy
#' @author Michael Mansfield
#'
#' @family Reducing functions
#'
#' @examples
#' exampleColinear5 |> wholeGenomeClassification()
#'
#' @export

wholeGenomeClassification <- function(gb, coa = coalesce_contigs(gb)) {
  # Isolated regions are found identical in the gb and coa objects.
  isol        <- granges(gb[gb %in% coa])
  if(length(isol) != 0)
    isol$type <- "isolated alignment"

  # Collinear regions are unique to the coa objets.
  col         <- granges(gb[ ! gb %in% coa])
  if(length(col) != 0)
    col$type  <- "collinear alignment"

  # Regions not represented in gb are unaligned.
  unal        <- cleanGaps(gb)

  # Unaligned regions are 'mapped' if they are flanked by collinear regions.
  # Unmapped regions are not represented in coa objects.
  unmap       <- cleanGaps(coa)

  # Breakpoint regions are unaligned and unmapped.
  brk         <- granges(unal[unal %in% unmap])
  if(length(brk) != 0)
    brk$type  <- "breakpoint region"

  # Bridge regions are unaligned, but they are mapped.
  bri        <- granges(unal[ ! unal %in% unmap ])
  if(length(bri) != 0)
    bri$type   <- "bridge region"

  wholeGenome <- c(isol, brk, col, bri) |> sort(ignore.strand = TRUE)
  wholeGenome$type <- factor(wholeGenome$type,
                             levels = c("isolated alignment", "collinear alignment",
                                        "breakpoint region", "bridge region"))
  wholeGenome
}
