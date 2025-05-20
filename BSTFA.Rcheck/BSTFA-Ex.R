pkgname <- "BSTFA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BSTFA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BSTFA")
### * BSTFA

flush(stderr()); flush(stdout())

### Name: BSTFA
### Title: Reduced BSTFA function
### Aliases: BSTFA

### ** Examples

data(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords)



cleanEx()
nameEx("BSTFAfull")
### * BSTFAfull

flush(stderr()); flush(stdout())

### Name: BSTFAfull
### Title: Full BSTFA function
### Aliases: BSTFAfull

### ** Examples

data(utahDataList)
out <- BSTFA.full(ymat=TemperatureVals, dates=Dates, coords=Coords)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
