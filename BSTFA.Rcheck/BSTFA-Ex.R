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
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=1000)



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
out <- BSTFAfull(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=50)



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
nameEx("compute_summary")
### * compute_summary

flush(stderr()); flush(stdout())

### Name: compute_summary
### Title: Print computation summary
### Aliases: compute_summary

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
compute_summary(out)



cleanEx()
nameEx("convergence_diag")
### * convergence_diag

flush(stderr()); flush(stdout())

### Name: convergence_diag
### Title: Check effective sample size and geweke diagnostic
### Aliases: convergence_diag

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
convergence_diag(out)



cleanEx()
nameEx("map_spatial_param")
### * map_spatial_param

flush(stderr()); flush(stdout())

### Name: map_spatial_param
### Title: Plot a map of interpolated spatially-dependent parameter values.
### Aliases: map_spatial_param

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
map_spatial_param(out, parameter="slope", map=TRUE, state=TRUE, location='utah', fine=50)



cleanEx()
nameEx("plot_annual")
### * plot_annual

flush(stderr()); flush(stdout())

### Name: plot_annual
### Title: Plot annual curve
### Aliases: plot_annual

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot_annual(out, location=1)



cleanEx()
nameEx("plot_factor")
### * plot_factor

flush(stderr()); flush(stdout())

### Name: plot_factor
### Title: Plot the factors
### Aliases: plot_factor

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot_factor(out, factor=1:4, together=TRUE)



cleanEx()
nameEx("plot_location")
### * plot_location

flush(stderr()); flush(stdout())

### Name: plot_location
### Title: Plot a location's time series of estimated/interpolated values.
### Aliases: plot_location

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot_location(out, location=1, pred.int=FALSE)



cleanEx()
nameEx("plot_spatial_param")
### * plot_spatial_param

flush(stderr()); flush(stdout())

### Name: plot_spatial_param
### Title: Plot the spatially-dependent parameter for in-sample locations.
### Aliases: plot_spatial_param

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot_spatial_param(out, parameter="slope")



cleanEx()
nameEx("plot_trace")
### * plot_trace

flush(stderr()); flush(stdout())

### Name: plot_trace
### Title: Plot trace plots
### Aliases: plot_trace

### ** Examples

data(utahDataList)
attach(utahDataList)
out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
plot_trace(out, parameter="beta", param.range=1)



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
