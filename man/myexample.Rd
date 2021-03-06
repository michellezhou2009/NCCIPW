% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myexample.R
\docType{data}
\name{myexample}
\alias{myexample}
\title{Example NCC data}
\format{
A list with  elements, which are:
\describe{
\item{data}{a \code{data.frame} with the following columns: \code{time}, \code{status}, \code{marker1}, and \code{marker2}}
\item{id}{ID of the subjects}
\item{case}{binary values indicating whether each subject is a sub-cohort case}
\item{control}{binary values indicating whether the subject is a sub-cohort control}
\item{m0}{the number of controls selected for each case}
\item{matchid}{a matrix in which the number of rows is equal to the number of selected sub-cohort cases, and the number of columns is \code{m0}+1. In each row, the first element is the ID of a case and the remaining elements are the IDs of the controls selected for the given case}
\item{Mdat}{a matrix of two matching variables}
\item{aM}{a 2-dimensional vector of matching contraints}
}
}
\usage{
myexample
}
\description{
An example of NCC data including the survival outcomes,  \code{time} (censored event time) and \code{status} (censoring indicator variable) and two markers: \code{marker1} and \code{marker2}.
}
\keyword{datasets}
