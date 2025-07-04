% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flagInversions.R
\name{leftInversionGaps}
\alias{leftInversionGaps}
\title{Isolate the left-side gaps in inversions}
\usage{
leftInversionGaps(gb)
}
\arguments{
\item{gb}{A \code{\link{GBreaks}} object.}
}
\value{
Returns a \code{\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}} object representing the left-side gaps in the
\code{GBreaks} object.
}
\description{
Isolate the left-side gaps in inversions
}
\examples{
leftInversionGaps(flagInversions(exampleInversion))

}
\seealso{
Other Inversion functions: 
\code{\link{filterDoubleInversions}()},
\code{\link{filterInversions}()},
\code{\link{flagDoubleInversions}()},
\code{\link{flagInversions}()},
\code{\link{flagPairs}()},
\code{\link{flipInversions}()},
\code{\link{showInversions}()}
}
\concept{Inversion functions}
