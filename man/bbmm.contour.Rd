\name{bbmm.contour}
\alias{bbmm.contour}
\title{bbmm.contour - calculate and plot density contours for a Brownian bridge 
    movement model. 
}

\description{Calculates and (optionally) plots density contours (e.g. 99%) for a 
    utilization distribution. This function was built for the BBMM package and 
    the output from \code{\link{brownian.bridge}}, but could be used to produce
    contours for any probabilty density grid. The contours are determined by 
    integration using the trapezoidal rule.
}

\usage{
bbmm.contour(x, levels, locations = NULL, plot = TRUE)
}

\arguments{
  \item{x}{An object of class bbmm, which is a list with vectors X, Y and probabilities.
}

 \item{levels}{A numeric vector of desired contour levels (e.g., 90, 95, 99).
}

  \item{locations}{(optional) data frame or matrix of x and y coordinates for original
  location data used to estimate the Brownian bridge and probability density surface.
}
  \item{plot}{Logical. If true, the countour plot is drawn.
}

}

\value{
An list with four components. Use names(fit), where the call was fit <- contour.bbmm(...), 
to see names of all returned components.

    Components of the returned object are as follows:
  \item{Contour}{Requested contour level(s).}
  \item{Z}{Z values that separate the specified contour level(s).}

}

\author{
Ryan Nielson, Hall Sawyer, and Trent McDonald (www.west-inc.com)

Maintainer: Ryan Nielson <rnielson@west-inc.com>
}

\seealso{
 \code{\link{brownian.bridge}}
}

\examples{
data(locations)
BBMM <- brownian.bridge(x=locations$x, y=locations$y, 
            time.lag=locations$time.lag, location.error=20, 
            cell.size=50)
bbmm.summary(BBMM)
contours <- bbmm.contour(BBMM, levels=c(seq(50, 90, by=10), 95, 99), 
        locations=locations, plot=TRUE)
print(contours)

# Create data.frame indicating cells within the 50% contour and export as Ascii Grid
contour.50 <- data.frame(x = BBMM$x, y = BBMM$y, probability = BBMM$probability) 
contour.50 <- contour.50[contour.50$probability >= contours$Z[1],] 
contour.50$in.out <- 1

# Not run: require(maptools) 
# Not run: m <- SpatialPixelsDataFrame(points = contour.50[c("x", "y")], data=contour.50) 
# Not run: m <- as(m, "SpatialGridDataFrame") 
# Not run: writeAsciiGrid(m, "contour.asc", attr=ncol(m)) 

}

\keyword{ ~contour }
\keyword{ ~Ascii }
