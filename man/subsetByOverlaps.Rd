% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenomicBreaks.R
\name{subsetByOverlaps_GBreaks}
\alias{subsetByOverlaps_GBreaks}
\alias{subsetByOverlaps,GBreaks,GBreaks-method}
\title{subsetByOverlaps method for \code{GBreaks} objects}
\usage{
subsetByOverlaps_GBreaks(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...
)

\S4method{subsetByOverlaps}{GBreaks,GBreaks}(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...
)
}
\arguments{
\item{x}{A \code{GBreaks} object.}

\item{ranges}{Another \code{GBreaks} object.}

\item{maxgap}{tbd}

\item{minoverlap}{tbd}

\item{type}{tbd}

\item{invert}{tbd}

\item{...}{etc}
}
\value{
tbd
}
\description{
This is a \code{subsetByOverlaps} method for \code{\link{GBreaks}} objects, that will run
\code{\link[GenomicRanges:findOverlaps-methods]{GenomicRanges::subsetByOverlaps}} on its \emph{target} and \emph{query} ranges and
will return a new \code{GBreaks} object.
}
\note{
\code{range} unconditionally ignores strand in \code{GBreaks} objects.
}
\examples{
subsetByOverlaps(exampleColinear3, exampleColinear3)
}
\seealso{
Other Bioconductor API functions: 
\code{\link{GBreaks-class}},
\code{\link{getSeq}()},
\code{\link{pairwiseAlignment}()},
\code{\link{range_GBreaks}()},
\code{\link{reverse}()}
}
\concept{Bioconductor API functions}
