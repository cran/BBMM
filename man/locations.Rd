\name{locations}
\alias{locations}
\docType{data}
\title{Example data for the Brownian bridge movement model function in package BBMM.}
\description{Seventeen GPS locations from a female mule deer. Locations were recorded every 2 hours.}
\usage{data(locations)}
\format{
  A data frame with 17 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{A numeric vector of X coordinates}
    \item{\code{y}}{A numeric vectorof Y coordinates}
    \item{\code{time.lag}}{A numeric vector of time differences between each successive observation.}
  }
}

\examples{
data(locations)
BBMM = brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag, location.error=20, 
            cell.size=50)
bbmm.summary(BBMM)
contours = bbmm.contour(BBMM, levels=seq(50, 90, by=10), 
        locations=locations, plot=TRUE)
print(contours)

}
\keyword{datasets}
