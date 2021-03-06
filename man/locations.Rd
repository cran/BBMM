\name{locations}
\alias{locations}
\docType{data}
\title{
Example data for the Brownian bridge movement model function in package BBMM.
}
\description{
25 GPS locations from a female mule deer. Locations were recorded every 2 hours.
}
\usage{data(locations)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{A numeric vector of x coordinates (UTM)}
    \item{\code{y}}{A numeric vectorof y coordinates (UTM)}
    \item{\code{time.lag}}{A numeric vector of time differences (min) between each successive observation.}
  }
}
\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag[-1], location.error=20, 
            cell.size=50)
bbmm.summary(BBMM)
}
\keyword{datasets}
