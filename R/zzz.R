
.onLoad <- function(lib, pkg) {
    vg190 = compareVersion(paste(version$major, 
                                 version$minor, sep="."), "1.9.0")
    
    if (vg190 < 0) {
        if(!require(ctest))
            warning("Could not load package ctest")
    } else {
        if(!require(stats))                 
            warning("Could not load package stats")      
    }
}
