bbmm.summary <-
function(x){
#
#   print method for patch occupancy objects
#
    cat("\nBrownian motion variance : ", x[[1]], fill=TRUE) 
    cat("Size of Grid : ", length(x$x), "Cells", fill=TRUE)
    cat( "Grid Cell Size : ", abs(x$x[1]-x$x[2]), fill=TRUE)

}

