\name{bbmm.contour}
\alias{bbmm.contour}
\title{
Create density contours from a Brownian bridge movement model
}
\description{
Calculates and (optionally) plots density contours from a utilization distribution. 
}
\usage{
bbmm.contour(x, levels, locations = NULL, plot = TRUE)
}
\arguments{
	\item{x}{An object of class bbmm, which is a list with vectors \code{x}, \code{y} and \code{probabilities}.
	}

	\item{levels}{A numeric vector of desired contour levels (e.g., \code{c(95, 99)}).
	}

	\item{locations}{(optional) data frame or matrix of \code{x} and \code{y} coordinates for original
	  location data used to estimate the Brownian bridge and utilization distribution surface.
	}
	
	\item{plot}{Logical. If \code{TRUE}, the countour plot is drawn.
	}
}
\value{
An list with four components. 

    Components of the returned object are as follows:
  \item{contour}{Requested contour level(s).}
  \item{z}{values that separate the specified contour level(s).}
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
contours <- bbmm.contour(BBMM, levels=c(95, 99), locations=locations, plot=TRUE)

# Create a shapefile with contour lines
# Not run: library(maptools)
# Not run: library(raster)
# Not run: out <- data.frame(x=BBMM$x,y=BBMM$y,z=BBMM$probability)
# Make sure the data is properly projected
# Not run: out.raster <- rasterFromXYZ(out,crs=CRS("+proj=utm +zone=12 +datum=WGS84"),digits=2) 
# Not run: raster.contour <- rasterToContour(out.raster,levels=contours$Z) 
# Not run: raster.contour <- spChFIDs(raster.contour,paste(c(95, 99),"% Contour Line",sep="")) 
# Not run: library(rgdal) # can't be loaded before CRS call in line above
# Not run: writeOGR(obj=raster.contour,dsn=".",layer="BBMM",driver="ESRI Shapefile")

}

\keyword{ ~contour }
\keyword{ ~Ascii }
