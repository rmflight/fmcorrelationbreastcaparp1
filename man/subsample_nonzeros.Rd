% Generated by roxygen2 (4.0.2): do not edit by hand
\name{subsample_nonzeros}
\alias{subsample_nonzeros}
\title{subsample non-zeros}
\usage{
subsample_nonzeros(data, data_columns, log_transform = TRUE,
  non_zero = "either", n_points = 1000)
}
\arguments{
\item{data}{a \code{DataFrame}}

\item{data_columns}{which columns to use}

\item{log_transform}{whether or not to log-transform the data}

\item{non_zero}{"either" one of the columns, or "both"}

\item{n_points}{how many points to sample}
}
\value{
data.frame
}
\description{
Takes a \code{DataFrame} instance and for the two columns indicated, returns a subset of points that are non-zero in one or both of
the variables.
}

