% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flagAll.R
\name{flagAll}
\alias{flagAll}
\title{Compute all flags}
\usage{
flagAll(gb)
}
\arguments{
\item{gb}{A \code{GBreaks} object.}
}
\value{
Returns the \code{GBreaks} object with an extra \code{flag} metadata column.
}
\description{
Flag ranges that are at the beginning of a colinear duplet, or an inversion
or a translocation triplet.
}
\examples{
flagAll(exampleInversion)

}
\seealso{
Other Flagging functions: 
\code{\link{flagColinearAlignments}()},
\code{\link{flagDoubleInversions}()},
\code{\link{flagInversions}()},
\code{\link{flagLongShort}()},
\code{\link{flagPairs}()},
\code{\link{flagTranslocations}()}
}
\concept{Flagging functions}
