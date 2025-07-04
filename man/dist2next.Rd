% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist2next.R
\name{dist2next}
\alias{dist2next}
\alias{dist2next,GRanges-method}
\alias{dist2next,GBreaks-method}
\title{Distances to next ranges}
\usage{
dist2next(x, step = 1, ignore.strand = FALSE)

\S4method{dist2next}{GRanges}(x, step = 1, ignore.strand = FALSE)

\S4method{dist2next}{GBreaks}(x, step = 1, ignore.strand = FALSE)
}
\arguments{
\item{x}{A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}} or a \code{\link{GBreaks}} object.}

\item{step}{Distance to the \ifelse{html}{\out{n<sup>th</sup>}}{\eqn{2n^th}}
block (default: first).}

\item{ignore.strand}{Calculate distance for ranges on different strands.}
}
\value{
For \code{GRanges}, returns the object with one extra metadata colums,
\code{dist}, and for \code{GBreaks} two extra columns \code{tdist} and \code{qdist}, containing
the distance to the next range in the reference and query genomes
respectively.
}
\description{
Calculates the distance to the next range in the \emph{target} and \emph{query} genome.
}
\details{
The distance is defined by the \code{\link[GenomicRanges:nearest-methods]{GenomicRanges::distance}} function.  Note
that because the last range has not next neighbor, the last value is set
to \code{Inf} arbitrarily.
}
\examples{
dist2next(exampleInversion)
dist2next(granges(exampleInversion))
dist2next(exampleInversion, ignore.strand = TRUE)
dist2next(exampleInversion - 20, ignore.strand = TRUE)
dist2next(exampleInversion, 2)

}
\seealso{
Other Colinearity functions: 
\code{\link{GOC}()},
\code{\link{bridgeRegions}()},
\code{\link{chain_contigs}()},
\code{\link{coalesce_contigs}()},
\code{\link{filterColinearRegions}()},
\code{\link{flagColinearAlignments}()},
\code{\link{flagPairs}()}
}
\concept{Colinearity functions}
