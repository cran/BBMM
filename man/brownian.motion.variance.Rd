\name{brownian.motion.variance}
\alias{brownian.motion.variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Brownian motion variance}
\description{Estimates Brownian motion variance for a Brownian bridge movement model.}

\usage{
brownian.motion.variance(n.locs, time.lag, location.error, x, y)
}
\arguments{
  \item{n.locs}{Number of observed locations in time and space.
}
  \item{time.lag}{Vector of time differences (usually in minutes) between 
    successive locations. The length(time.lag) = length(x) - 1. 
}
  \item{location.error}{The standard deviation of normally distributed location errors (single value or vector of 
  1 value for each observation).}
  \item{x}{Vector of x coordinates (meters) of locations, ordered in time.}
  \item{y}{Vector of y coordinates (meters) of locations, ordered in time.}

}

\details{Users will generally not call this function directly. It is sourced by 
    \code{\link{brownian.bridge}}. A maximum likelihood estimate is obtained
    using \code{\link{optimize}}.
}

\value{
The estimate of Brownian motion variance, or diffusion coefficient related to 
    the mobility of the object, defined as deqn{sigma_{m}^2} in Horne et al. (2007)).
}

\references{
Horne, J. S., E. O. Garton, S. M. Krone, and J. S. Lewis. 2007. 
    Analyzing animal movements using Brownian bridges. 
    Ecology 88:2354-2363.

Sawyer, H., M. J. Kauffman, R. M. Nielson, and J. S. Horne. 2009. 
    Identifying and prioritizing ungulate migration routes for 
    landscpe-level conservation. Ecological Applications 19:2016-2025.

}

\author{
Ryan Nielson, Hall Sawyer, and Trent McDonald (www.west-inc.com)

Maintainer: Ryan Nielson <rnielson@west-inc.com>

}

\seealso{\code{\link{brownian.bridge}}
}

\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag, location.error=20, 
            cell.size=50)

}

\keyword{ ~Brownian }
\keyword{ ~motion }
\keyword{ ~variance }
