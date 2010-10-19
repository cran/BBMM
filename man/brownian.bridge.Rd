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
brownian.bridge(x, y, time.lag, location.error, area.grid = NULL, cell.size = NULL, time.step = 10)

}

\arguments{
  \item{x}{Vector of x coordinates (in meters) of locations, ordered in time.}
  \item{y}{Vector of y coordinates (in meters) of locations, ordered in time.}
  \item{time.lag}{Vector of time differences (in minutes) between 
  successive locations (length(time.lag) = length(x) - 1).}
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
  discretizing time into arbitrarily small intervals of time.step. The default 
  is 10 units (same as time.lag). Longer time.step speeds up estimation, 
  but reduces precision.}

}

\details{
This is the main routine for estimating a Brownian bridge. It calls 
brownian.motion.variance to estimate the Brownian motion variance via 
maximum likelihood, and then a Fortran routine to calculate the 
probabilities of use across the area.grid.

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
            time.lag=locations$time.lag, location.error=20, 
            cell.size=50)

}

\keyword{ ~Brownian }
\keyword{ ~bridge }
\keyword{ ~movement }
\keyword{ ~GPS }
\keyword{ ~UD }
