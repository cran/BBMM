.First.lib<-function(libname, pkgname){
    library.dynam("BBMM", pkgname)
    cat(paste( "Brownian brige movement model", "\n\nWritten by: \n\tRyan Nielson (rnielson@west-inc.com), \n\tHall Sawyer, \n\tand Trent McDonald\n\t(www.west-inc.com)\n") )
}
