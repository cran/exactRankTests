
.onLoad <- function(lib, pkg) {
    packageStartupMessage(paste(" Package", sQuote("exactRankTests"), 
        "is no longer under development.\n",
        "Please consider using package", sQuote("coin"), "instead.\n"))
}
