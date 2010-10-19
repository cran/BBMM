\name{BBMM}
\alias{BBMM}
\docType{package}
\title{BBMM - Brownian bridge movement model}
\description{
This package fits a Brownian bridge movement model to observed locations 
    in space and time. This is a continuous-time stochastic model of movement 
    in which the probability of being in an area during the time of observation 
    is conditioned on starting and ending locations. A BBMM is typically fit 
    to animal location data obtained by a Global Positioning System (GPS) 
    or Very High Frequency (VHF) device. The model provides an empirical 
    estimate of the movement path of an animal using discrete location 
    data obtained at relatively short time intervals. This package allows the user
    to specify the grid cell size and extent over which the probability of use
    will be estimated.

}

\details{
\tabular{ll}{
Package: \tab BBMM\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2010-10-19\cr
License: \tab GNU General Public License\cr
LazyLoad: \tab yes\cr

}

List of routines:

\code{\link{brownian.bridge}}            Work-horse routine for estimating a Brownian bridge movement model.

\code{\link{brownian.motion.variance}}   Estimates the Brownian motion variance.

\code{\link{bbmm.summary}}               Summary method for bbmm objects.

\code{\link{bbmm.contour}}               Contour plot for bbmm objects.
 
}

\author{
Ryan Nielson, Hall Sawyer, and Trent McDonald  (www.west-inc.com)

Maintainer: Ryan Nielson <rnielson@west-inc.com>

}

\references{
Horne, J. S., E. O. Garton, S. M. Krone, and J. S. Lewis. 2007. 
    Analyzing animal movements using Brownian bridges. 
    Ecology 88:2354-2363.

Sawyer, H., M. J. Kauffman, R. M. Nielson, and J. S. Horne. 2009. 
    Identifying and prioritizing ungulate migration routes for 
    landscpe-level conservation. Ecological Applications 19:2016-2025.

}

\keyword{ package }
