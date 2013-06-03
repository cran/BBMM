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
  \item{time.lag}{Vector of time differences (in minutes) between successive locations 
  where \code{length(time.lag)} is the same as \code{length(x)-1}. 
  \code{time.lag[1]} is the length of time between locations 1 and 2, and 
  \code{time.lag[2]} is the difference between locations 2 and 3, and so on.}
  \item{location.error}{The standard deviation of normally distributed location errors (single value or vector of 
  1 value for each observation).}
  \item{area.grid}{(optional) Matrix or data frame of x and y coordinates of cell 
  center ponts on a rectangular grid that defines the area in which to 
  estimate probability of use. If missing, a grid is created by expanding
  the range of x and y by 1 standard deviation and using cell.size.}
  \item{cell.size}{(optional) Cell size for \code{area.grid}, if \code{area.grid} not 
  provided. Must specify either \code{area.grid} or \code{cell.size}.}
  \item{time.step}{(optional) The Brownian bridge probability density function 
  must be integrated to find the fraction of time spent in each region. While the 
  probability density function cannot be integrated, it can be approximated by 
  discretizing time into arbitrarily small intervals according to time.step. The default 
  is 10 units (same as \code{time.lag}). A longer \code{time.step} speeds up estimation, 
  but reduces precision.}
  \item{max.lag}{(optional) Maximum time lag (same units as \code{time.lag}) between successful locations 
  to use in calculating the Brownian motion variance and probability of use. This can be important 
  if some scheduled locations are missing, either via an unsuccessful fix attempt to locate the individual/object or 
  because of a lack of effort. Including large time gaps in sequence of locations can artificially 
  inflate/deflate the Brownian motion variance and potentially bias estimates of probability of 
  use (to what degree is unknown). If two successive locations are \code{>max.lag}, then a Brownian bridge
  is not estimated between those two locations. The default sets the \code{max.lag} equal to the maximum time 
  difference in \code{time.lag+1}.}
}
\details{
This is the main routine for estimating a Brownian bridge. It calls 
brownian.motion.variance to estimate the Brownian motion variance via 
maximum likelihood and then calculates the probabilities of use across the \code{area.grid}. Larger data sets
and larger grids require more computing time, which can be a few of hours on a 32-bit PC or just a 
fraction of an hour on a 64-bit PC running R x64.  
}
\value{
An object (list) of class("bbmm") with four components. 

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
Ryan Nielson, Hall Sawyer, and Trent McDonald (WEST, Inc., \url{www.west-inc.com})

Maintainer: Ryan Nielson \email{rnielson@west-inc.com}
}
\seealso{\code{\link{brownian.motion.variance}}, \code{\link{bbmm.summary}}, 
  \code{\link{bbmm.contour}} 
}
\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag[-1], location.error=20, 
            cell.size=50)

# To export an Ascii grid with probabilities (UD).
# First delete any x and y coords that have probability of use < 0.00000001.
# 	This saves some computing time and reduces the size of the output Ascii file.
# Not run: x <- BBMM$x[BBMM$probability >= 0.00000001]
# Not run: y <- BBMM$y[BBMM$probability >= 0.00000001]
# Not run: z <- BBMM$probability[BBMM$probability >= 0.00000001]
# Not run: tmp <- data.frame(x, y, z) 
# Not run: library(maptools)
# Not run: m <- SpatialPixelsDataFrame(points = tmp[c("x", "y")], data=tmp)
# Not run: m <- as(m, "SpatialGridDataFrame")
# Not run: write.asciigrid(m, "BBMM.asc", attr=3)

}
\keyword{ ~Brownian }
\keyword{ ~bridge }
\keyword{ ~movement }
\keyword{ ~GPS }
\keyword{ ~UD }
