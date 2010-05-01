brownian.bridge <-
function(x, y, time.lag, location.error, area.grid=NULL, cell.size=NULL, time.step=10){

    # Input: 
    #   x, y = vectors of coordinates of ordered animal locations, in UTMs
    #   time.lag = vector of time differences (in minutes) between successive 
    #           locations. length(time.lag) = length(x) - 1. 
    #   location.error = either a mean location error or a vector of errors 
    #                   (1 for each location)
    #   area.grid = matrix or data frame of X and Y coords for brownian bridge. 
    #           If missing, defaults to minimum/maximum x/y minus/plus 10% of the 
    #           range of x/y. If area.grid is provided, cell sizes must be square and uniform. 
    #   cell.size = cell size for grid, if grid not provided
    #   time.step = 
    # Output: 
    #   UD = list with estimated Brownian.Motion.Variance, X, Y, and Z, 
    #       where X and Y are grid center point coordinates and Z is the estimated
    #       probability of use with sum(Z) = 1.0.
    
    if(is.null(x) | is.null(y) | (length(x) != length(y))) {
        stop("data is missing or unequal number of x and y coordinates")
    }
    if(is.null(location.error)) stop("must specify 'location.error'")
    if(is.null(area.grid) & is.null(cell.size)) {
        stop("'area.grid' or 'cell.size' must be specified")
    }
    if(!is.null(area.grid) & is.null(cell.size)){
        cell.size = abs(area.grid[1,1] - area.grid[2,1]) 
    }
    if(is.null(area.grid) & !is.null(cell.size)){
        range.x = range(x)
        range.y = range(y)
        min.grid.x = round(range.x[1] - 1*sd(x))
        max.grid.x = round(range.x[2] + 1*sd(x))
        min.grid.y = round(range.y[1] - 1*sd(y))
        max.grid.y = round(range.y[2] + 1*sd(y))
        X=seq(min.grid.x, max.grid.x, cell.size) 
        Y=seq(min.grid.y, max.grid.y, cell.size)
        area.grid = merge(X, Y)
    }
    
    if(length(location.error) == 1){
        location.error = rep(location.error, length(x))
    }
    
    n.locs = length(x)

    BMvar = brownian.motion.variance(n.locs, time.lag, location.error, x, y)
    BMvar = rep(BMvar, times=length(x))

    # Use 10 units (generally minutes) as default.
    if(is.null(time.step)) time.step = 10

    grid.size = nrow(area.grid)
    probability = rep(0, grid.size)
    T.Total = sum(time.lag)

    bbmm = vector("list", 4)
    names(bbmm) = c("Brownian motion variance", "x", "y", "probability")
    class(bbmm) = "bbmm"
    
    #dyn.load("BBMM.dll")

    ans <- .Fortran("BBMM",
                    as.integer(n.locs), 
                    as.integer(grid.size), 
                    as.double(time.lag), 
                    as.double(T.Total), 
                    as.double(x), 
                    as.double(y), 
                    as.double(BMvar), 
                    as.double(location.error), 
                    as.double(area.grid$x), 
                    as.double(area.grid$y), 
                    as.double(time.step), 
                    as.double(probability),
                    PACKAGE="BBMM")
    
    bbmm[[4]] = ans[[12]]    

    #dyn.unload("BBMM.dll")

    bbmm[[1]] = BMvar[1]
    bbmm[[2]] = area.grid$x
    bbmm[[3]] = area.grid$y

    return(bbmm)

}

