
.onLoad <- function(lib, pkg) {
    if(!require(ctest))
        warning("Could not load package ctest")
}
