#' Convert `GRanges` object to `dna_seg` format
#'
#' The [`genoPlotR::dna_seg`] class represents _DNA segments_ using
#' a `data.frame` format.
#'
#' @param gr A [`GRanges`] object.
#'
#' @param ... Extra arguments passed to `dna_seg()`.
#'
#' @returns Returns a `genoPlotR::dna_seq` or `NULL` if the input object has a
#' length of zero.
#'
#' @author Charles Plessy
#'
#' @importFrom genoPlotR dna_seg
#'
#' @examples
#' inv       <- GRanges(c("XSR:101-180:+", "XSR:201-300:-",  "XSR:320-400:+"))
#' inv$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' gr2dna_seg(inv)
#' gr2dna_seg(inv$query)
#'
#' @family genoPlotR functions
#' @family converter functions
#'
#' @export

gr2dna_seg <- function (gr, ...) {
  # DNA segments need unique names
  if (length(gr) == 0) return(NULL)
  if (is.null(names(gr))) names(gr) <- as.character(gr)
  # Unstranded ranges are not recognised
  strand(gr[strand(gr) == "*"]) <- "+"
  df <- data.frame(
    name   = names(gr),
    start  = start(gr),
    end    = end(gr),
    strand = strand(gr)
  )
  genoPlotR::dna_seg(df, gene_type = "blocks", ...)
}

#' Convert `GBreaks` object to a list of `dna_seg` objects
#'
#' @returns Returns a `list` of two [`genoPlotR::dna_seg`] objects, respectively
#' for the _target_ ranges the _query_ ranges.
#'
#' @param gb A [`GRanges`] or a [`GBreaks`] object.
#'
#' @param ... Extra arguments passed to `dna_seg()`.
#'
#' @author Charles Plessy
#'
#' @importFrom genoPlotR dna_seg
#'
#' @examples
#' inv       <- GRanges(c("XSR:101-180:+", "XSR:201-300:-",  "XSR:320-400:+"))
#' inv$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' inv <- GBreaks(inv)
#' gb2dna_seg(inv)
#'
#' @family genoPlotR functions
#' @family converter functions
#'
#' @export

gb2dna_seg <- function(gb, ...) {
  target <- gr2dna_seg(gb, ...)
  query  <- gr2dna_seg(gb$query, ...)
  # By convention in GenomicBreaks, strand information refers to the query
  s <- target$strand
  target$strand <- query$strand
  query$strand <- s
  list(target=target, query=query)
}

#' Convert `GBreaks` object to a `comparison` object
#'
#' @returns Returns a [`genoPlotR::comparison`] object.
#'
#' @note Note that the `comparison` objects assume that the ranges are on the
#' same sequence feature.  It is your responsibility that the input `GBreaks`
#' object also conforms to this assumption.
#'
#' @param gb A [`GBreaks`] object.
#' @param color (optional) A vector of colors, of same length as `gb`.
#' @param ignore.strand Ignore strand information?
#'
#' @author Charles Plessy
#'
#' @importFrom genoPlotR as.comparison
#'
#' @examples
#' inv       <- GRanges(c("XSR:101-180:+", "XSR:201-300:-",  "XSR:320-400:+"))
#' inv$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' inv <- GBreaks(inv)
#' gb2comp(inv)
#'
#' @family genoPlotR functions
#' @family converter functions
#'
#' @export

gb2comp <- function(gb, color = NULL, ignore.strand = FALSE) {
  dna_segs <- gb2dna_seg(gb)
  df <- data.frame(
    start1 = dna_segs$target$start,
    end1   = dna_segs$target$end,
    start2 = dna_segs$query$start,
    end2   = dna_segs$query$end
  )
  cmp <- genoPlotR::as.comparison(df)
  cmp$direction <- 1L
  # If color information is absent, the genoPlotR::plot_gene_map function will
  # bug and fail to plot comparisons on direction `-1`.
  cmp$color <- "darksalmon"
  if(isFALSE(ignore.strand)) {
    onMinus <- decode(strand(gb) == "-")
    cmp$direction[onMinus] <- -1
    cmp$start2[onMinus]    <- dna_segs$query$end  [onMinus]
    cmp$end2  [onMinus]    <- dna_segs$query$start[onMinus]
    cmp$col   [onMinus]    <- "cornflowerblue"
  }
  if(!is.null(color)) {
    cmp$col <- color
  }
  cmp
}

#' Plots a pair of chromosomes from a `GBreaks` object.
#'
#' Plots the mapping between a pair of chromosomes using the
#' [`genoPlotR::plot_gene_map`] function.  One sequence feature (usually a
#' chromosome) is selected from the _target_ genome.  The corresponding sequence
#' feature on the _query_ genome is either chosen automatically (being the one
#' with the largest fraction of the mappings) or given as a parameter.
#'
#' @note In this kind of plot, the ranges from the selected sequence on the
#' _target_ genome that are not mapped to the _query_ sequence are not
#' displayed.
#'
#' @returns Plots to the active device and returns and `oma::oma_layout` object
#' like the [`genoPlotR::plot_gene_map`] function.
#'
#' @param gb A [`GBreaks`] object.
#' @param chrT A sequence name on the _target_ genome.  Defaults to the first
#'        _sequence level_ of the `gb` object.
#' @param chrQ (Optional) A sequence name on the _query_ genome.  Defaults to
#'        the longest cumulative match on `chrT`.
#' @param ... Further arguments are passed to `plot_gene_map`.
#'
#' @author Charles Plessy
#'
#' @importFrom genoPlotR plot_gene_map
#' @importFrom BiocGenerics `%in%`
#'
#' @examples
#' plotApairOfChrs(exampleInversion)
#'
#' @family genoPlotR functions
#' @family plot functions
#' @family Structural variants
#'
#' @export

plotApairOfChrs <- function(gb, chrT=NULL, chrQ=NULL, ...) {
  # If needed, guess name of target seqlevel
  if(is.null(chrT))
    chrT <- seqlevelsInUse(gb)[1]

  # Object restricted to the target seqlevel
  gb.chrT <- gb[seqnames(gb) == chrT]

  # If needed, guess name of query seqlevel
  if(is.null(chrQ))
    chrQ <- tapply(width(gb.chrT$query), seqnames(gb.chrT$query), sum) |> sort() |> tail(1) |> names()

  # Object restricted to the query seqlevel
  gb.chrQ <- gb[seqnames(gb$query) == chrQ]

  # Object restricted to both seqlevels
  roi <- gb.chrT[seqnames(gb.chrT$query) == chrQ]

  # Ranges of gb.chrT not in the roi
  extraT <- granges(gb.chrT[! gb.chrT %in% roi])

  # Ranges of gb.chrQ not in the roi
  extraQ <- gb.chrQ$query[! gb.chrQ$query %in% roi$query]

  dsList        <- gb2dna_seg(roi)
  # Append extra blocks with `rbind` to handle `NULL` objects.
  dsList$target <- rbind(dsList$target, gr2dna_seg(extraT, col = "lightblue", fill = "lightblue"))
  dsList$query  <- rbind(dsList$query,  gr2dna_seg(extraQ, col = "lightblue", fill = "lightblue"))

  compList <- list(gb2comp(roi))

  genoPlotR::plot_gene_map(dsList, compList, ...)
}
