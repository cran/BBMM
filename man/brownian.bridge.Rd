\name{brownian.bridge}
\alias{brownian.bridge}
\alias{Brownian}
\alias{bridge}
\alias{movement}
\alias{model}
\title{Brownian bridge movement model}
\description{Estimate a Brownian bridge model of movement in which the 
probability of a mobile object being in an area is conditioned on starting 
and ending locations. The model provides an empirical estimate of a movement 
path using discrete location data obtained at relatively short time 
intervals. The Brownian bridge probability density connecting each pair of 
successive locations is an estimate of the relative time spent in an area 
during the time interval between those locations. 

}
\usage{
brownian.bridge(x, y, time.lag, location.error, area.grid = NULL, 
	cell.size = NULL, time.step = 10, max.lag = NULL)
}

\arguments{
  \item{x}{Vector of x coordinates (in meters) of locations, ordered in time.}
  \item{y}{Vector of y coordinates (in meters) of locations, ordered in time.}
  \item{time.lag}{Vector of time differences (in minutes) between 
  successive locations (length(time.lag) = length(x) - 1). I.e., time.lag[1] is the length of time between
  locations 1 and 2, and time.lag[2] is the difference between locations 2 and 3, and so on.}
  \item{location.error}{The standard deviation of normally distributed location errors (single value or vector of 
  1 value for each observation).}
  \item{area.grid}{(optional) Matrix or data frame of x and y coordinates of cell 
  center ponts on a rectangular grid that defines the area in which to 
  estimate probability of use. If missing, a grid is created by expanding 
  the range of x and y by 1 standard deviation and using cell.size.}
  \item{cell.size}{(optional) Cell size for area.grid, if area.grid not 
  provided. Must specify either area.grid or cell.size.}
  \item{time.step}{(optional) The Brownian bridge probability density function 
  must be integrated to find the fraction of time spent in each region. While the 
  probability density function cannot be integrated, it can be approximated by 
  discretizing time into arbitrarily small intervals according to time.step. The default 
  is 10 units (same as time.lag). A longer time.step speeds up estimation, 
  but reduces precision.}
  \item{max.lag}{(optional) Maximum time lag (same units as time.lag) between successful locations 
  to use in calculating the brownian motion variance and probability of use. This can be important 
  if some data are missing, either via an unsuccessful attempt to locate the individual/object or 
  because of a lack of effor. Including large time gaps in sequence of locations can artificially 
  inflate/deflate the Brownian motion variance and potentially bias estimates of probability of 
  use (to what degree is unknown). If two successive locations are > max.lag, then a Brownian bridge
  is not estimated between those two locations. The default sets the max.lag equal to the maximum time 
  difference in time.lag + 1.}
}

\details{
This is the main routine for estimating a Brownian bridge. It calls 
brownian.motion.variance to estimate the Brownian motion variance via 
maximum likelihood and then calculates the probabilities of use across the area.grid. Larger data sets
and larger grids require more computing time, which can be a few of hours on a 32-bit PC or just a 
fraction of an hour on a 64-bit PC.  

}

\value{
An object (list) of class("bbmm") with four components. Use bbmm.summary to 
print a summary. Use names(fit), where the call was fit <- brownian.bridge(...), 
to see names of all returned components.

    Components of the returned object are as follows:
  \item{Brownian motion variance}{Estimated Brownian motion variance.}
  \item{x}{Vector of x coordinates for grid cells.}
  \item{y}{Vector of y coordinates for grid cells.}
  \item{probability}{Estimated proportion of time spent in each grid cell.}

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

\seealso{\code{\link{brownian.motion.variance}}, \code{\link{bbmm.summary}}, 
  \code{\link{bbmm.contour}} 

}
  
\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag[-1], location.error=20, 
            cell.size=50)

}

\keyword{ ~Brownian }
\keyword{ ~bridge }
\keyword{ ~movement }
\keyword{ ~GPS }
\keyword{ ~UD }