pkgname <- "BSTFA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BSTFA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BSTFA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BSTFA")
### * BSTFA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BSTFA
### Title: Reduced BSTFA function
### Aliases: BSTFA

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BSTFA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("BSTFAfull")
### * BSTFAfull

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BSTFAfull
### Title: Full BSTFA function
### Aliases: BSTFAfull

### ** Examples

## Not run: 
##D #Example below not run; even the ten iterations will take a minute or two to run.
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFAfull(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=10)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BSTFAfull", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeLogLik")
### * computeLogLik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeLogLik
### Title: Compute log-likelihood
### Aliases: computeLogLik

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D loglik <- computeLogLik(out, addthin=2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeLogLik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_summary")
### * compute_summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_summary
### Title: Print computation summary
### Aliases: compute_summary

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D compute_summary(out)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("convergence_diag")
### * convergence_diag

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: convergence_diag
### Title: Check effective sample size and geweke diagnostic
### Aliases: convergence_diag

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D convergence_diag(out)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("convergence_diag", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("map_spatial_param")
### * map_spatial_param

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: map_spatial_param
### Title: Plot a map of interpolated spatially-dependent parameter values.
### Aliases: map_spatial_param

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D map_spatial_param(out, parameter='slope', map=TRUE, state=TRUE, location='utah', fine=50)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("map_spatial_param", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_annual")
### * plot_annual

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_annual
### Title: Plot annual/seasonal behavior at a specific location.
### Aliases: plot_annual

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D plot_annual(out, location=1)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_annual", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_factor")
### * plot_factor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_factor
### Title: Plot the temporally-dependent factors.
### Aliases: plot_factor

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D plot_factor(out, factor=1:4, together=TRUE)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_factor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_location")
### * plot_location

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_location
### Title: Plot a location's time series of estimated/predicted values.
### Aliases: plot_location

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D plot_location(out, location=1, pred.int=FALSE)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_location", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_spatial_param")
### * plot_spatial_param

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_spatial_param
### Title: Plot the spatially-dependent parameter for in-sample locations.
### Aliases: plot_spatial_param

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D plot_spatial_param(out, parameter='slope')
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_spatial_param", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_trace")
### * plot_trace

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_trace
### Title: Plot trace plots
### Aliases: plot_trace

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D plot_trace(out, parameter='beta', param.range=1)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_trace", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predictBSTFA")
### * predictBSTFA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predictBSTFA
### Title: Estimate/predict values of the time series at a specific
###   location.
### Aliases: predictBSTFA

### ** Examples

## Not run: 
##D data(utahDataList)
##D attach(utahDataList)
##D out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
##D loc1means <- predictBSTFA(out, location=1, pred.int=FALSE)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predictBSTFA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
