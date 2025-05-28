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
attach(utahDataList)
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
attach(utahDataList)
out <- BSTFA.full(ymat=TemperatureVals, dates=Dates, coords=Coords)



cleanEx()
nameEx("check.convergence")
### * check.convergence

flush(stderr()); flush(stdout())

### Name: check.convergence
### Title: Check effective sample size and geweke diagnostic
### Aliases: check.convergence

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
check.convergence(out)



cleanEx()
nameEx("computation.summary")
### * computation.summary

flush(stderr()); flush(stdout())

### Name: computation.summary
### Title: Print computation summary
### Aliases: computation.summary

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
computation.summary(out)



cleanEx()
nameEx("computeLogLik")
### * computeLogLik

flush(stderr()); flush(stdout())

### Name: computeLogLik
### Title: Compute log-likelihood
### Aliases: computeLogLik

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
loglik <- computeLogLik(out, addthin=2)



cleanEx()
nameEx("plot.annual")
### * plot.annual

flush(stderr()); flush(stdout())

### Name: plot.annual
### Title: Plot annual curve
### Aliases: plot.annual

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.annual(out, location=1)



cleanEx()
nameEx("plot.factor")
### * plot.factor

flush(stderr()); flush(stdout())

### Name: plot.factor
### Title: Plot the factors
### Aliases: plot.factor

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.factor(out, factor=1:4, together=T)



cleanEx()
nameEx("plot.grid")
### * plot.grid

flush(stderr()); flush(stdout())

### Name: plot.grid
### Title: Plot the spatially-dependent parameter for in-sample locations.
### Aliases: plot.grid

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.grid(out, parameter="slope")



cleanEx()
nameEx("plot.location")
### * plot.location

flush(stderr()); flush(stdout())

### Name: plot.location
### Title: Plot a location's time series of estimated/interpolated values.
### Aliases: plot.location

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.location(out, location=1, pred.int=FALSE)



cleanEx()
nameEx("plot.map")
### * plot.map

flush(stderr()); flush(stdout())

### Name: plot.map
### Title: Plot a map of interpolated spatially-dependent parameter values.
### Aliases: plot.map

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.map(out, parameter="slope", map=T, state=T, location='utah', fine=50)



cleanEx()
nameEx("plot.trace")
### * plot.trace

flush(stderr()); flush(stdout())

### Name: plot.trace
### Title: Plot trace plots
### Aliases: plot.trace

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot.trace(out, parameter="beta", param.range=1)



cleanEx()
nameEx("predictBSTFA")
### * predictBSTFA

flush(stderr()); flush(stdout())

### Name: predictBSTFA
### Title: Prediction
### Aliases: predictBSTFA

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
loc1means <- predictBSTFA(out, location=1, pred.int=FALSE)



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
