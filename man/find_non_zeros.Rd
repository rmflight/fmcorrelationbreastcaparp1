% Generated by roxygen2 (4.0.2): do not edit by hand
\name{find_non_zeros}
\alias{find_non_zeros}
\title{find non-zeros}
\usage{
find_non_zeros(data, data_columns, log_transform = TRUE,
  non_zero = "either")
}
\arguments{
\item{data}{the data we are working with}

\item{data_columns}{which columns to use}

\item{log_transform}{do a log transformation on the data before identifying zeros}

\item{non_zero}{"either" or "both", how much data *must* be non-zero}
}
\value{
indices into the original data that are not zero
}
\description{
given a \code{data.frame}, find the non-zero entries in each case
}

