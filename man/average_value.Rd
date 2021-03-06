% Generated by roxygen2 (4.0.2): do not edit by hand
\name{average_value}
\alias{average_value}
\title{weighted average by width}
\usage{
average_value(granges, use_column = "mcols.signal")
}
\arguments{
\item{granges}{the \code{GRanges} object}

\item{use_column}{which column from \code{mcols} has the signal to average}
}
\value{
numeric
}
\description{
given a \code{GRanges} object and an \code{mcols} column, a weighted average of the \code{mcols} column will
be calculated where the weight will be the fractional number of bases indicated by the width of each entry in
the \code{GRanges} object.
}

