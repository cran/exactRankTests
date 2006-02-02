
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

    cat("  Package", sQuote("exactRankTests"), "is no longer under development.\n",
        " Please consider using package", sQuote("coin"), "instead.\n")
}
