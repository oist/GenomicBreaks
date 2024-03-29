% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flagColinearAlignments.R
\name{filterColinearRegions}
\alias{filterColinearRegions}
\title{Filter colinear regions}
\usage{
filterColinearRegions(gb, rename = TRUE)
}
\arguments{
\item{gb}{A \code{\link{GBreaks}} object processed with \code{\link{flagColinearAlignments}}.}

\item{rename}{Replace range names by their numeric order before subsetting.}
}
\value{
Returns the \code{GBreaks} object in which all ranges that are not
the central part of an inversion triplet have been discarded.  If the object
was missing the \code{colinear} metadata column, return the object after
discarding all of its ranges.
}
\description{
Filter colinear regions
}
\examples{
filterColinearRegions(flagColinearAlignments(exampleColinear))

}
\seealso{
Other Colinearity functions: 
\code{\link{GOC}()},
\code{\link{bridgeRegions}()},
\code{\link{chain_contigs}()},
\code{\link{coalesce_contigs}()},
\code{\link{dist2next}()},
\code{\link{flagColinearAlignments}()},
\code{\link{flagPairs}()}
}
\concept{Colinearity functions}
