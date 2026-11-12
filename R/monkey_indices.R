#' Monkey indices
#'
#' Twelve summary statistics of a pairwise genome alignment that were not designed
#' to represent a biological fact.  They were retained from a screen of 5,420
#' functions assembled mechanically from the columns of the alignment table
#' and from generic descriptive statistics ("monkey typing"), for the following
#' properties on 1,910 alignments of 16 clades: reproducible between the
#' alignment of _A_ to _B_ and the alignment of _B_ to _A_, more alike between
#' genome pairs that share an ancestor than between random pairs, different
#' from the value obtained after shuffling the order and the strand of the
#' alignments, and not strongly correlated with the statistics that already
#' existed in this package, nor with each other (absolute rank correlation
#' below 0.6).  Their biological meaning, if any, is to be found.  Each function
#' is named after a famous monkey.
#'
#' @details
#'
#' All the functions take the alignments as loaded and sort them on the _target_
#' genome.  `mizaru_index()` then coalesces them with [`coalesce_contigs`] and
#' its default arguments; the others use them as they are.  They use the
#' following quantities.
#'
#'  - A _junction_ is a pair of consecutive alignments on the same _target_
#'    sequence.  A _same-query junction_ is a junction whose two alignments are
#'    also on the same _query_ sequence.
#'  - The _target gap_ of a junction is the number of bases between its two
#'    alignments on the _target_ genome, plus one; the _query gap_ of a
#'    same-query junction is the same on the _query_ genome, whatever the order
#'    and the strand of the two alignments.
#'  - The _gap contrast_ of a same-query junction is \eqn{(t - q) / (t + q)}
#'    where \eqn{t} and \eqn{q} are its target gap and its query gap.  It is
#'    positive when the two alignments are further apart on the _target_.
#'  - The _rank step_ of a same-query junction is the number of positions
#'    between its two alignments in the order of the alignments on the _query_
#'    genome (1 for neighbours), and its _direction_ is +1 when the two
#'    alignments are on the same strand and follow each other on the _query_
#'    as that strand implies, -1 when they follow each other in the opposite
#'    direction, and 0 when they are on opposite strands.
#'  - The _edge distance_ of an alignment on a genome is the distance from its
#'    centre to the nearest end of the sequence, divided by half the length of
#'    the sequence: 0 at the ends, 1 in the middle.  The _relative position_ is
#'    the position of the centre from 0 (first base) to 1 (last base).
#'  - _Sizes_ (widths, gaps, rank steps) are used on the logarithmic scale.
#'
#' The functions are the following.  The monkey is introduced first.
#'
#' \describe{
#'   \item{`wukong_index()`}{_Sun Wukong: the Monkey King of the Chinese novel
#'     Journey to the West._  Tail weight of the query gaps: on the logarithmic
#'     scale, range between the 10th and the 90th percentiles divided by the
#'     interquartile range.}
#'   \item{`abu_index()`}{_Abu: kleptomaniac pet monkey of the street thief in
#'     the Aladdin films._  Quartile skewness of the gap contrasts:
#'     \eqn{(Q_3 + Q_1 - 2 Q_2) / (Q_3 - Q_1)}.}
#'   \item{`albert_index()`}{_Albert II: rhesus macaque, first primate to reach
#'     space, in 1949._  Excess kurtosis of the gap contrasts.}
#'   \item{`marcel_index()`}{_Marcel: white-headed capuchin kept by Ross in the
#'     sitcom Friends._  Standard deviation of the absolute difference between
#'     the relative positions of an alignment on the two genomes.  The _query_
#'     position is counted from the other end of the sequence when most
#'     aligned bases between the two sequences are on the minus strand.}
#'   \item{`imo_index()`}{_Imo: Japanese macaque of Koshima who started washing
#'     sweet potatoes._  Correlation between the edge distance of an alignment
#'     on the _target_ genome and the edge distance of its neighbours on the
#'     _query_ genome.  Neighbours are the alignments before and after it on
#'     the same _target_ sequence; the two members of a junction are
#'     exchangeable.}
#'   \item{`rafiki_index()`}{_Rafiki: mandrill shaman who presents the newborn
#'     cub in The Lion King._  Rank correlation, over same-query junctions,
#'     between the query gap and the mean edge distance of the two alignments
#'     on the _query_ genome.}
#'   \item{`hanuman_index()`}{_Hanuman: monkey-faced Hindu deity, devoted
#'     companion of Rama in the Ramayana._  Rank correlation, over junctions,
#'     between the mean width of the two alignments and the absolute difference
#'     of their edge distances on the _query_ genome.}
#'   \item{`jack_index()`}{_Jack: chacma baboon who operated railway signals in
#'     South Africa in the 1880s._  Rank correlation, over same-query junctions,
#'     between the rank step and the direction.}
#'   \item{`able_index()`}{_Able: rhesus macaque, Miss Baker's companion on the
#'     1959 spaceflight._  Correlation, over same-query junctions, between the
#'     rank step and the absolute difference of the edge distances of the two
#'     alignments on the _query_ genome.}
#'   \item{`mizaru_index()`}{_Mizaru: the see-no-evil monkey carved at the
#'     Toshogu shrine in Nikko._  Rank correlation, over the junctions of the coalesced
#'     alignments, between the target gap and the absolute difference of the
#'     edge distances of the two coalesced alignments on the _query_ genome.}
#'   \item{`george_index()`}{_Curious George: inquisitive little monkey of the
#'     children's books by the Reys._  Correlation, over junctions, between the
#'     target gap and the mean score per aligned column of the two alignments.
#'     Stops with an error without the `score` and `aLength` metadata
#'     columns.}
#'   \item{`baker_index()`}{_Miss Baker: squirrel monkey who survived a
#'     suborbital spaceflight in 1959._  Correlation, weighted by the width of
#'     the alignments on the _target_ genome, between the number of gap
#'     openings per aligned column of an alignment and the geometric mean of
#'     the target gaps of the junctions that flank it.  Stops with an error
#'     without the `gapOpens` and `aLength` metadata columns.}
#' }
#'
#' Ranks are computed within each vector before restricting it to the
#' junctions used: in `rafiki_index()`, the mean edge distances are ranked
#' among all junctions.  Variances are population variances.  Quantiles are
#' those of [`stats::quantile`] with default settings.
#'
#' The results do not depend on the order of the sequences nor on their
#' orientation: reversing a sequence with [`reverse`] does not change them.
#' They are not symmetric: the value on `swap(gb)` differs.  Note that
#' [`swap`] drops the metadata columns that `george_index()` and
#' `baker_index()` need.
#'
#' The functions that use edge distances or relative positions need sequence
#' lengths.  Missing lengths are replaced by the end of the last alignment of
#' each sequence, with [`forceSeqLengths`].
#'
#' @param gb A [`GBreaks`] object, not coalesced, in which all ranges are on
#'        the `+` or `-` strand.
#'
#' @returns A numeric value, or `NA` when it can not be computed: fewer than
#' three values, or no variation.
#'
#' @family Monkey indices
#'
#' @seealso The vignette `vignette("MonkeyIndices", package = "GenomicBreaks")`
#' reports how the functions were selected and how they behave on the
#' alignments of the screen.  The code of the screen is in
#' `system.file("scripts", "monkey_screen", package = "GenomicBreaks")`.
#'
#' @author Claude (Anthropic), from a screen requested by Charles Plessy.
#'
#' @examples
#' # Alignment of two yeast genomes
#' gb <- load_genomic_breaks(
#'   system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"))
#' wukong_index(gb)
#' abu_index(gb)
#' albert_index(gb)
#' marcel_index(gb)
#' imo_index(gb)
#' rafiki_index(gb)
#' hanuman_index(gb)
#' jack_index(gb)
#' able_index(gb)
#' mizaru_index(gb)
#'
#' # The orientation of the sequences does not matter.
#' marcel_index(reverse(forceSeqLengths(gb)))
#'
#' # The GFF3 format does not record the alignment length nor the number of gap
#' # openings, which the blasttab+ format does.  Lower bounds are used here
#' # for illustration.
#' gb$aLength  <- pmax(width(gb), width(gb$query))
#' gb$gapOpens <- as.integer(width(gb) != width(gb$query))
#' george_index(gb)
#' baker_index(gb)
#'
#' @name monkey_indices
#' @aliases monkey_indices
NULL

# Internal: the input must be a GBreaks object with oriented ranges.
.mk_check <- function(gb) {
  if (!inherits(gb, "GBreaks")) stop("`gb` must be a GBreaks object.")
  if (any(decode(strand(gb)) == "*"))
    stop("Monkey indices need the orientation of all ranges, ",
         "but some ranges are on the `*` strand.")
  invisible(gb)
}

# Internal: plain vectors describing the alignments sorted on the target, the
# junctions, and (if lengths = TRUE) positions relative to sequence ends.
.mk_blocks <- function(gb, lengths = FALSE) {
  .mk_check(gb)
  gb <- sort(gb, ignore.strand = TRUE)
  n  <- length(gb)
  x  <- list(
    n = n, gb = gb,
    tseq = as.integer(seqnames(gb)), qseq = as.integer(seqnames(gb$query)),
    ts = as.numeric(start(gb)),       te = as.numeric(end(gb)),
    qs = as.numeric(start(gb$query)), qe = as.numeric(end(gb$query)),
    s  = ifelse(decode(strand(gb)) == "-", -1, 1))
  x$tw <- x$te - x$ts + 1
  x$qw <- x$qe - x$qs + 1
  # Rank of each alignment in the order of the query genome
  x$qrank <- integer(n)
  x$qrank[order(x$qseq, x$qs, x$qe)] <- seq_len(n)
  # Junctions: index of their first alignment
  i <- if (n > 1) which(x$tseq[-n] == x$tseq[-1]) else integer(0)
  x$J   <- i
  x$sameq <- x$qseq[i] == x$qseq[i + 1]
  x$JQ  <- i[x$sameq]
  if (isTRUE(lengths)) {
    gb <- forceSeqLengths(gb)
    tL <- as.numeric(seqlengths(gb))[x$tseq]
    qL <- as.numeric(seqlengths(gb$query))[x$qseq]
    x$tpos  <- ((x$ts + x$te) / 2 - 0.5) / tL
    x$qpos  <- ((x$qs + x$qe) / 2 - 0.5) / qL
    x$tedge <- round(1 - abs(2 * x$tpos - 1), 12)
    x$qedge <- round(1 - abs(2 * x$qpos - 1), 12)
  }
  x
}

# Internal: gaps of the junctions (number of bases between the alignments + 1)
.mk_tgap  <- function(x, i = x$J)
  pmax(x$ts[i + 1] - x$te[i] - 1, 0) + 1
.mk_qstep <- function(x, k = x$JQ)
  pmax(pmax(x$qs[k], x$qs[k + 1]) - pmin(x$qe[k], x$qe[k + 1]) - 1, 0) + 1
.mk_gdiff <- function(x) {
  tg <- .mk_tgap(x, x$JQ)
  qg <- .mk_qstep(x)
  round((tg - qg) / (tg + qg), 12)
}

# Internal: TRUE when a sum of squares of n values of mean m is only
# floating-point noise.
.mk_zero <- function(ss, n, m) ss <= 1e-18 * n * (m * m + 1)

# Internal: Pearson correlation, optionally weighted; NA if fewer than three
# values or no variation.
.mk_cor <- function(a, b, w = NULL) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 3) return(NA_real_)
  w  <- if (is.null(w)) rep(1, length(a)) else w[ok]
  sw <- sum(w)
  ma <- sum(w * a) / sw
  mb <- sum(w * b) / sw
  sa <- sum(w * (a - ma)^2)
  sb <- sum(w * (b - mb)^2)
  if (.mk_zero(sa, sw, ma) || .mk_zero(sb, sw, mb)) return(NA_real_)
  sum(w * (a - ma) * (b - mb)) / sqrt(sa * sb)
}

# Internal: metadata column, with an error if it is missing.
.mk_col <- function(gb, name) {
  col <- if (name == "score") score(gb) else mcols(gb)[[name]]
  if (is.null(col))
    stop("This index needs the `", name, "` metadata column, which is missing.")
  col
}

#' @rdname monkey_indices
#' @export

wukong_index <- function(gb) {
  x <- .mk_blocks(gb)
  y <- log(.mk_qstep(x))
  if (length(y) < 3) return(NA_real_)
  q <- unname(quantile(y, c(.1, .25, .75, .9)))
  if (q[3] <= q[2]) return(NA_real_)
  (q[4] - q[1]) / (q[3] - q[2])
}

#' @rdname monkey_indices
#' @export

abu_index <- function(gb) {
  x <- .mk_blocks(gb)
  y <- .mk_gdiff(x)
  if (length(y) < 3) return(NA_real_)
  q <- unname(quantile(y, c(.25, .5, .75)))
  if (q[3] <= q[1]) return(NA_real_)
  (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1])
}

#' @rdname monkey_indices
#' @export

albert_index <- function(gb) {
  x <- .mk_blocks(gb)
  y <- .mk_gdiff(x)
  if (length(y) < 3) return(NA_real_)
  d <- y - mean(y)
  v <- mean(d^2)
  if (.mk_zero(v * length(y), length(y), mean(y))) return(NA_real_)
  mean(d^4) / v^2 - 3
}

#' @rdname monkey_indices
#' @export

marcel_index <- function(gb) {
  x <- .mk_blocks(gb, lengths = TRUE)
  if (x$n < 3) return(NA_real_)
  # Majority strand of each pair of sequences, by aligned target bases.  Exact
  # ties are broken with other properties that do not depend on the
  # orientation of the sequences.
  g    <- as.integer(factor(paste(x$tseq, x$qseq)))
  gmaj <- rep(1, max(g))
  todo <- rep(TRUE, max(g))
  for (w in list(x$tw, x$qw, x$tedge, x$qedge)) {
    if (!any(todo)) break
    S <- as.vector(rowsum(w * x$s, g, reorder = TRUE))
    gmaj[todo & S < 0] <- -1
    todo <- todo & S == 0
  }
  qpos <- ifelse(gmaj[g] > 0, x$qpos, 1 - x$qpos)
  y <- round(abs(x$tpos - qpos), 12)
  v <- mean((y - mean(y))^2)
  if (.mk_zero(v * length(y), length(y), mean(y))) return(0)
  sqrt(v)
}

#' @rdname monkey_indices
#' @export

imo_index <- function(gb) {
  x <- .mk_blocks(gb, lengths = TRUE)
  i <- x$J
  if (length(i) < 3) return(NA_real_)
  # The two members of a junction are exchangeable: one mean and one variance
  # for both, so that the reading direction does not matter.
  mt <- (sum(x$tedge[i]) + sum(x$tedge[i + 1])) / (2 * length(i))
  mq <- (sum(x$qedge[i]) + sum(x$qedge[i + 1])) / (2 * length(i))
  te <- x$tedge - mt
  qe <- x$qedge - mq
  st <- (sum(te[i]^2) + sum(te[i + 1]^2)) / 2
  sq <- (sum(qe[i]^2) + sum(qe[i + 1]^2)) / 2
  if (.mk_zero(st, length(i), mt) || .mk_zero(sq, length(i), mq)) return(NA_real_)
  (sum(te[i] * qe[i + 1]) + sum(te[i + 1] * qe[i])) / 2 / sqrt(st * sq)
}

#' @rdname monkey_indices
#' @export

rafiki_index <- function(gb) {
  x <- .mk_blocks(gb, lengths = TRUE)
  i <- x$J
  if (length(x$JQ) < 3) return(NA_real_)
  # Mean edge distance on the query, ranked among all junctions
  edge <- rank(round((x$qedge[i] + x$qedge[i + 1]) / 2, 12))[x$sameq]
  .mk_cor(rank(.mk_qstep(x)), edge)
}

#' @rdname monkey_indices
#' @export

hanuman_index <- function(gb) {
  x <- .mk_blocks(gb, lengths = TRUE)
  i <- x$J
  if (length(i) < 3) return(NA_real_)
  lw <- log(x$tw)
  .mk_cor(rank(round((lw[i] + lw[i + 1]) / 2, 12)),
          rank(round(abs(x$qedge[i] - x$qedge[i + 1]), 12)))
}

#' @rdname monkey_indices
#' @export

jack_index <- function(gb) {
  x <- .mk_blocks(gb)
  k <- x$JQ
  if (length(k) < 3) return(NA_real_)
  step <- x$qrank[k + 1] - x$qrank[k]
  .mk_cor(rank(abs(step)), rank(sign(step) * (x$s[k] + x$s[k + 1]) / 2))
}

#' @rdname monkey_indices
#' @export

able_index <- function(gb) {
  x <- .mk_blocks(gb, lengths = TRUE)
  k <- x$JQ
  if (length(k) < 3) return(NA_real_)
  .mk_cor(log(abs(x$qrank[k + 1] - x$qrank[k])),
          round(abs(x$qedge[k] - x$qedge[k + 1]), 12))
}

#' @rdname monkey_indices
#' @export

mizaru_index <- function(gb) {
  .mk_check(gb)
  # Sequence lengths are settled before coalescing, so that edge distances
  # refer to the same lengths as in the other indices.
  x <- .mk_blocks(coalesce_contigs(forceSeqLengths(gb)), lengths = TRUE)
  i <- x$J
  if (length(i) < 3) return(NA_real_)
  .mk_cor(rank(.mk_tgap(x)),
          rank(round(abs(x$qedge[i] - x$qedge[i + 1]), 12)))
}

#' @rdname monkey_indices
#' @export

george_index <- function(gb) {
  x <- .mk_blocks(gb)
  score   <- .mk_col(x$gb, "score")
  aLength <- .mk_col(x$gb, "aLength")
  i <- x$J
  if (length(i) < 3) return(NA_real_)
  sdens <- round(score / aLength, 12)
  .mk_cor(log(.mk_tgap(x)), round((sdens[i] + sdens[i + 1]) / 2, 12))
}

#' @rdname monkey_indices
#' @export

baker_index <- function(gb) {
  x <- .mk_blocks(gb)
  gapOpens <- .mk_col(x$gb, "gapOpens")
  aLength  <- .mk_col(x$gb, "aLength")
  i <- x$J
  # Mean, on the logarithmic scale, of the target gaps that flank each alignment
  lg  <- log(.mk_tgap(x))
  sum <- cnt <- numeric(x$n)
  sum[i]     <- sum[i] + lg
  cnt[i]     <- cnt[i] + 1
  sum[i + 1] <- sum[i + 1] + lg
  cnt[i + 1] <- cnt[i + 1] + 1
  bgap <- ifelse(cnt > 0, sum / cnt, NA_real_)
  .mk_cor(round(gapOpens / aLength, 12), bgap, w = x$tw)
}
