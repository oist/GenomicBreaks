% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeSeqLevels.R
\name{mergeSeqLevels}
\alias{mergeSeqLevels}
\title{Merge seqlevels in a larger one}
\usage{
mergeSeqLevels(gr, seqs, name)
}
\arguments{
\item{gr}{A \code{GRanges} object.}

\item{seqs}{A character vector of \code{\link{seqlevels}} from \code{gr}}

\item{name}{The name of the new sequence level to be added}
}
\value{
Returns a modified \code{GRanges} object in which the sequences have been
merged.  Its \code{\link{seqinfo}} has a new entry for the new level, and the old
levels are not removed.  If no \code{seqlengths} were present in the original
object, they are arbitrarily set as the maximal end value for each \code{seqlevel}.
}
\description{
For scaffolding or plotting purposes, it may be useful to merge some
sequences into larger ones.
}
\note{
Be careful that in some cases it is needed to "flip" the sequence
feature with \code{\link{reverse}} before merging, for instance when colinearity is
with its reverse strand.
}
\examples{
gb       <- GRanges(c("XSR:101-180:+", "XSR:201-300:+",  "XSR:320-400:+"))
gb$query <- GRanges(c( "S1:101-200",      "S2:1-100",    "S3:1-100"))
seqlengths(gb$query) <- c(200, 100, 100)
genome(gb$query) <- "GenomeX"
isCircular(gb$query) <- rep(FALSE, 3)
seqinfo(gb$query)
gb <- GBreaks(gb)
gb$query <- mergeSeqLevels(gb$query, c("S2", "S3"), "Scaf1")
gb
seqinfo(gb$query)

mergeSeqLevels(gb, seqlevelsInUse(gb), "AllMerged")

}
\seealso{
Other modifier functions: 
\code{\link{bridgeRegions}()},
\code{\link{coalesce_contigs}()},
\code{\link{forceSeqLengths}()},
\code{\link{guessSeqLengths}()},
\code{\link{reverse}()},
\code{\link{swap}()}
}
\author{
Charles Plessy
}
\concept{modifier functions}
