\name{bbmm.summary}
\alias{bbmm.summary}
\title{
Summary of a Brownian bridge movement model fit
}
\description{
A function used to produce a result summary from a Brownian bridge movement model.
}
\usage{
bbmm.summary(x)
}
\arguments{
  \item{x}{
	An object of class 'bbmm', which is a result of \code{brownian.bridge(...)}.
   }
}
\value{
The estimated Brownian motion variance, along with the number and size of 
grid cells over which the Brownian bridge was estimated.
}
\author{
Ryan Nielson, Hall Sawyer, and Trent McDonald (WEST, Inc., \url{www.west-inc.com})

Maintainer: Ryan Nielson \email{rnielson@west-inc.com}
}
\seealso{
\code{\link{brownian.bridge}}
}
\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag[-1], location.error=20, 
            cell.size=50)
bbmm.summary(BBMM)
}
\keyword{ ~brownian.bridge }
