.First.lib <- function(lib, pkg) {
  require(ctest)
  library.dynam("exactRankTests", pkg, lib)
}
