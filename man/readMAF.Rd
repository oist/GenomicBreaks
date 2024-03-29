% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{readMAF}
\alias{readMAF}
\title{Read a MAF file}
\usage{
readMAF(inputFileName)
}
\arguments{
\item{inputFileName}{The name of the file to read}
}
\value{
a \code{list} object with containg coordinates of the alignments in both
genomes and other information such as alignment width and number of matches.
}
\description{
Reads a pairwise genome alignment in MAF format.  The file can be plain
text or compressed with \code{gzip}.
}
\details{
Known limitations: Does not expand shell metacharacters.  Trusts blindly
file extension to determine compression.  Does not perform any validation on
the file format.  Assumes that the score comes first in the 'a' lines.
}
