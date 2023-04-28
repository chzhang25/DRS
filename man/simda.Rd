\name{simda}
\alias{simda}
\docType{data}
\title{
A simulated data with five longitudinal risk factors and one survival outcome subject to both independent and dependent censoring
}
\description{
A simulated data with five longitudinal risk factors and one survival outcome subject to both independent and dependent censoring
}
\usage{data("simda")}
\format{
  A data frame with five longitudinal risk factors and one survival outcome subject to both independent and dependent censoring.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{w}}{a continuous covariate with fixed value for each subject, e.g., age}
    \item{\code{v.time}}{visit time at which the longitudinal risk factors are measured}
    \item{\code{obsT}}{observed event time}
    \item{\code{delta}}{observed event indicator: 0 = independent censor, 1 = event, 2 = dependent censor}
    \item{\code{Lt1 - Lt5}}{five continuous longitudinal risk factors }
  }
}
\examples{
data(simda)
}
\keyword{datasets}
