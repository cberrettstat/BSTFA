## ----setup1, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache=FALSE,
  dev = "ragg_png",
  dpi=200,
  fig.width = 8
)

## ----setup2, echo=FALSE, include=FALSE----------------------------------------
library(knitr)
devtools::load_all()
library(kableExtra)

## ----fig-utahTemps, echo=FALSE, fig.height=3.5, fig.width=6, fig.cap='Mean-centered 30-day average daily minimum temperatures for three Utah weather stations (Moab, green; Canyonlands, orangee; Logan, purple) from January 1999 through December 2001.', fig.align="center"----

mycols <- RColorBrewer::brewer.pal(8, 'Dark2')
colnames(utahDataList$TemperatureVals) = utahDataList$Locations
window=1050:1090
plot(y=utahDataList$TemperatureVals[window,which(utahDataList$Locations=='MOAB')],
     x=utahDataList$Dates[window], type = 'l',
     main = "Measurements at Three Locations",
     xlab = "",
     ylab = "",
     ylim=c(-18,18),
     cex.main=1, col=mycols[1], lwd=2)
lines(y=utahDataList$TemperatureVals[window,which(utahDataList$Locations=='CANYONLANDS.THE.NECK')],
     x=utahDataList$Dates[window], col=mycols[2], lwd=2)
lines(y=utahDataList$TemperatureVals[window,utahDataList$Locations=='LOGAN.UTAH.ST.UNIV'], 
     x=utahDataList$Dates[window], col=mycols[3], lwd=2)
legend("bottomleft", legend=c("Moab", "Canyonlands", "Logan"), col=mycols[1:3], lty=1, lwd=2, bg="white", cex=.9)

## ----tab-required, results='asis', label="tab-required", echo=FALSE-----------
kable(
  data.frame(
    "Argument" = c("\\texttt{ymat}", "\\texttt{dates}", "\\texttt{coords}"),
    "Description" = c(
      "A matrix of response values. Each row should represent a point in time; each column should represent a specific location. Missing values should be recorded as NA.",
      "A vector of dates of length \\texttt{nrow(ymat)}. The model functions will transform this vector into 'doy' (day of year) using \\texttt{lubridate::yday()}. Thus, this vector must either be a lubridate or string object with year-month-day format.", 
      "A matrix of coordinate values with number of rows equal to \\texttt{ncol(ymat)} and 2 columns; if using longitude/latitude, longitude should be the first column."
    ),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Required arguments for the \\texttt{BSTFA} and \\texttt{BSTFAfull} functions."
) %>%  
  column_spec(2, width = "4.5in")

## ----BSTFAFunction, echo=FALSE, include=FALSE, cache=FALSE--------------------
# This code runs the function and is used for the plotting functions, but isn't shown in the vignette
set.seed(23416)
bstfa.output = BSTFA(ymat=utahDataList$TemperatureVals,
                     dates=utahDataList$Dates,
                     coords=utahDataList$Coords,
                     factors.fixed = c(144,89,129,78), # Specific fixed factor locations for Utah data set
                     iters=50,
                     burn=10,
                     verbose=TRUE,
                     save.output=FALSE)

## ----reducedmod, eval=FALSE, cache=FALSE--------------------------------------
# bstfa.output = BSTFA(ymat=utahDataList$TemperatureVals,
#                      dates=utahDataList$Dates,
#                      coords=utahDataList$Coords,
#                      verbose=FALSE)

## ----fullmod, eval=FALSE------------------------------------------------------
# bstfa_full.output = BSTFAfull(ymat=utahDataList$TemperatureVals,
#                          dates=utahDataList$Dates,
#                          coords=utahDataList$Coords,
#                          verbose=FALSE)

## ----predall, eval=F----------------------------------------------------------
# loc = 1 # Alpine, Utah in our data set
# preds = predictBSTFA(bstfa.output,
#                      location = loc,
#                      type='all',
#                      pred.int=TRUE,
#                      ci.level=c(0.025,0.975))

## ----plotpredobs, eval=F------------------------------------------------------
# loc = 1 # Alpine, Utah in our data set
# preds = predictBSTFA(bstfa.output,
#                      location = loc,
#                      type='mean')

## ----plotpred, eval=F---------------------------------------------------------
# loc = data.frame('Longitude' = 111.41, 'Latitude' = 38.29) # Torrey, Utah
# preds_new = predictBSTFA(bstfa.output,
#                          location = loc,
#                          type='all',
#                          pred.int=TRUE,
#                          ci.level=c(0.025,0.975))

## ----fig-location, fig.cap="Posterior mean temperature values (black line) at an in-sample location, Alpine, Utah, with 95% posterior predictive bounds (gray bands), and observed temperature measurements (gray circles).", fig.align="center", fig.height=3.5, fig.width=6----
loc = 1 # Alpine, Utah in our data set
plot_location(bstfa.output,
              location=loc,
              type='mean',
              pred.int=TRUE,
              uncertainty=TRUE,
              ci.level=c(0.025,0.975),
              xrange=c('1959-01-01', '1979-01-01'),
              truth=TRUE)

## ----fig-prediction, fig.cap="Posterior mean temperature values (black line) at an out-of-sample location, Torrey, Utah, with 95% posterior predictive bounds (gray bands).", fig.align="center",fig.height=3.5, fig.width=6----
loc = data.frame('Longitude' = 111.41, 'Latitude' = 38.29) # Torrey, Utah
plot_location(bstfa.output,
              location=loc,
              type='mean',
              pred.int=TRUE,
              uncertainty=TRUE,
              ci.level=c(0.025,0.975),
              xrange=c('1959-01-01', '1979-01-01'),
              truth=FALSE)

## ----tab-plotting, results='asis', label="tab-plotting", echo=FALSE-----------
kable(
  data.frame(
    "Function" = c("\\texttt{plot\\_location}", "\\texttt{plot\\_annual}", "\\texttt{plot\\_spatial\\_param}", "\\texttt{map\\_spatial\\_param}", "\\texttt{plot\\_factor}"),
    "Description" = c(
      "Plot predicted response variable at a specific location (either observed or unobserved) for a specified time range. Credible or prediction interval bands for a given probability (default is 95\\%) can be included. Uses base R for plotting.",
      "Plot the estimated annual seasonal behavior at a specific location (either observed or unobserved). Credible interval bands for a given probability (default is 95\\%) can be included. Uses base R for plotting.", 
      "Plot the estimated spatially-dependent linear slope or specific factor loading (the user specifies the parameter of interest) at all observed locations. Credible interval bounds for a given probability (default is 95\\%) can also be plotted. Uses \\texttt{ggplot2} for plotting.", 
      "Plot the interpolated spatially-dependent linear slope or specific factor loading (the user specifies the parameter of interest) on a grid of unobserved locations. Credible interval bounds for a given probability (default is 95\\%) can also be plotted. Contains arguments to import and plot the grid on a map using functions from the \\texttt{sf} package. Uses \\texttt{ggplot2} for plotting.", 
      "Plot the estimated factors, either individually or all together. Credible interval bands for a given probability (default is 95\\%) can be included. Uses base R for plotting."
    ),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Plotting/visualization functions in the \\texttt{BSTFA} package.", 
  linesep = ""
) %>%  
  column_spec(2, width = "4.5in")

## ----fig-season, fig.cap="Posterior mean annual seasonal behavior (black line) with 95% credible interval (gray band), and observed data (open gray circles) plotted by day of year.", fig.align="center",fig.height=3.5, fig.width=5----
plot_annual(bstfa.output,
            location=1,
            years='one')

## ----fig-obsslope, fig.align="center", fig.cap="Posterior mean linear change in time (slope) of temperatures at observed locations.", fig.width=3.5, fig.height=3.5----
plot_spatial_param(bstfa.output,
          type='mean',
          parameter='slope',
          yearscale=TRUE)

## ----fig-obsloading1, fig.cap="Posterior mean estimates of loadings for factor 1 at observed locations.  The circled red dot shows the location of the fixed factor loading.", fig.align="center", fig.width=4, fig.height=3.5----
plot_spatial_param(bstfa.output,
          type='mean',
          parameter='loading',
          loadings=1)

## ----fig-factorexample2, fig.cap="Posterior mean (black line) and 95% credible interval (gray band) of the first factor across the observed time period.", fig.align="center", fig.height=3.5, fig.width=6----
plot_factor(bstfa.output,
            together=FALSE,
            include.legend=FALSE,
            factor=1,
            type='mean')

## ----fig-factorexample1, fig.cap="Posterior mean (solid lines) and 95% credible intervals (light colored bands) for the four estimated factors across the observed time period.", fig.align="center", fig.height=3.75, fig.width=6----
plot_factor(bstfa.output,
            together=TRUE,
            include.legend=TRUE,
            type='mean')

## ----fig-mapexample1, cache=FALSE, warning=FALSE, fig.cap="Posterior mean linear change in temperatures across time for locations across the observed space.", fig.align="center", fig.width=3.5----
map_spatial_param(bstfa.output,
         parameter='slope',
         yearscale=TRUE,
         type='mean',
         map=TRUE,
         state=TRUE,
         location='utah',
         fine=50)

## ----tab-computation, results='asis', label="tab-computation", echo=FALSE-----
kable(
  data.frame(
    "$n$" = seq(100, 500, by=100),
    "\\texttt{BSTFA}, 8 loading bases" = c(0.016, 0.029, 0.051, 0.082, 0.12), 
    "\\texttt{BSTFA}, 20 loading bases" = c(0.017, 0.031, 0.05, 0.078, 0.119), 
    "\\texttt{BSTFA}, 50 loading bases" = c(0.021, 0.036, 0.056, 0.086, 0.126), 
    "\\texttt{BSTFAfull}" = c(.481, 1.292, 1.693, 2.7, 4.186),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Computation time in seconds per MCMC iteration for simulated data with $n$ locations (first column) and $T=300$ time points fit with the \\texttt{BSTFA} function for different number of Fourier bases for the loadings (8, 20, and 50, indicated by the 2nd, 3rd, and 4th columns) all with $R_t = 60$ Fourier bases for the factors, compared to the \\texttt{BSTFAfull} function (last column).",
  linesep = ""
)

## ----fig-load3, fig.cap="Comparison of estimated loadings for factor 3 using BSTFAfull (left), BSTFA with 8 spatial bases (center), and BSTFA with 6 spatial bases (right).", out.width="6.5in", echo=F, fig.align="center"----
knitr::include_graphics("loading_comparison3.png")

## ----fig-factor2, fig.cap="Comparison of estimated factors for factor 2 using BSTFAfull (left), BSTFA with 200 temporal bases (center), and BSTFA with 126 bases (right).", out.width="6.5in", echo=F, fig.align="center"----
knitr::include_graphics("factor_comparison2.png")

## ----comptime-----------------------------------------------------------------
compute_summary(bstfa.output)

## ----tab-linear, results='asis', label="tab-linear", echo=FALSE---------------
kable(
  data.frame(
    "Argument" = c("\\texttt{linear}", "\\texttt{beta}", "\\texttt{alpha.prec}", "\\texttt{tau2.gamma}", "\\texttt{tau2.phi}"),
    "Default Value" = c("\\texttt{TRUE}", "\\texttt{NULL}", "\\texttt{1e-5}", "\\texttt{2}", "\\texttt{1e-6}"),
    "Description" = c(
      "TRUE/FALSE value indicating whether the linear component is included in the model.",
      "Vector of starting values for $\\boldsymbol{\\beta}$ of length $n \\times 1$; if none is supplied, realistic starting values are calculated.",
      "Value on the diagonal of the precision matrix $\\mathbf{A}$.",
      "Value of the shape parameter, $\\gamma$, for the prior of the variance, $\\tau^2_\\beta$.",
      "Value of the rate parameter, $\\phi$, for the prior of the variance, $\\tau^2_\\beta$."
    ),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Arguments to \\texttt{BSTFA} and \\texttt{BSTFAfull} associated with the linear component."
) %>%  
  column_spec(3, width = "4in")

## ----tab-seasonal, results='asis', label="tab-seasonal", echo=FALSE-----------
kable(
  data.frame(
    "Argument" = c("\\texttt{seasonal}", "\\texttt{xi}", "\\texttt{n.seasn.knots}"),
    "Default Value" = c("\\texttt{TRUE}", "\\texttt{NULL}", "\\texttt{7}"),
    "Description" = c(
      "TRUE/FALSE value indicating whether the seasonal component should be included in the model.",
      "Vector of starting values for $\\boldsymbol{\\xi}$ of length $un \\times 1$; if none is supplied, realistic starting values are calculated.",
      "Value representing the value of $u$, the number of circular B-spline knots."),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Arguments to \\texttt{BSTFA} and \\texttt{BSTFAfull} associated with the seasonal component."
) %>%  
  column_spec(3, width = "4in")

## ----tab-factorcomp, results='asis', label="tab-factorcomp", echo=FALSE-------
kable(
  data.frame(
    "Argument" = c("\\texttt{factors}", "\\texttt{Fmat}", "\\texttt{Lambda}", "\\texttt{factors.fixed}", "\\texttt{n.factors}", "\\texttt{plot.factors}"),
    "Default Value" = c("\\texttt{TRUE}", "\\texttt{NULL}", "\\texttt{NULL}", "\\texttt{NULL}", "\\texttt{min(4, ceiling(ncol(ymat)/20))}", "\\texttt{FALSE}"),
    "Description" = c(
      "TRUE/FALSE value indicating whether the factor analysis component should be included in the model.",
      "Matrix of starting values for $\\mathbf{F}$ of dimension $T \\times L$; if none is supplied, realistic starting values are calculated.",
      "Matrix of starting values for $\\boldsymbol{\\Lambda}$ of dimension $n \\times L$; if none is supplied, realistic starting values are calculated.", 
      "Vector of indices (representing specific columns of \\texttt{ymat}) indicating locations to fix for the factors. If no vector is supplied, fixed factor locations are optimally chosen according to distance and amount of non-missing data. If this vector is supplied, \\texttt{n.factors = length(factors.fixed)}.", 
      "Number of factors to fit. If the number of locations is greater than 80, the function will always fit 4 factors unless otherwise specified in the \\texttt{factors.fixed} argument.",
      "TRUE/FALSE value indicating whether to provide a base R plot of the fixed factor locations."),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Arguments to \\texttt{BSTFA} and \\texttt{BSTFAfull} associated with the factor analysis component.",
  linesep = ""
) %>%  
  column_spec(2, width = "1.7in") %>%
  column_spec(3, width="3in")

## ----tab-resid, results='asis', label="tab-resid", echo=FALSE-----------------
kable(
  data.frame(
    "Argument" = c("\\texttt{sig2}", "\\texttt{sig2.gamma}", "\\texttt{sig2.phi}"),
    "Default Value" = c("\\texttt{NULL}", "\\texttt{2}", "\\texttt{1e-5}"),
    "Description" = c(
      "A starting value for $\\sigma^2$; if none is supplied, the starting value will be the variance of the non-missing values of \\texttt{ymat}.",
      "Value of the shape parameter, $\\gamma_\\sigma$, for the prior of the residual variance, $\\sigma^2$.",
      "Value of the rate parameter, $\\phi_\\sigma$, for the prior of the overall residual variance, $\\sigma^2$."
    ),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Arguments to \\texttt{BSTFA} and \\texttt{BSTFAfull} associated with the residual."
) %>%  
  column_spec(3, width = "4in")

## ----fig-fourierplots, fig.width=6, fig.cap="Fourier bases for the two-dimensional space.  Top row represents the sine bases (r = 1, 3, 5; left, center, right, respectively).  Bottom row represents the corresponding cosine bases.", fig.align='center'----
plot_fourier_bases(utahDataList$Coords,
                   R=6,
                   plot.3d=TRUE,
                   freq.lon = 4*diff(range(utahDataList$Coords[,1])),
                   freq.lat = 4*diff(range(utahDataList$Coords[,2])))

## ----tab-bases, results='asis', label="tab-bases", echo=FALSE-----------------
kable(
  data.frame(
    "Argument" = c("\\texttt{spatial.style}", "\\texttt{n.spatial.bases}", "\\texttt{load.style}", "\\texttt{n.load.bases}", "\\texttt{freq.lon}", "\\texttt{freq.lat}", "\\texttt{n.temp.bases}", "\\texttt{freq.temp}", "\\texttt{knot.levels}", "\\texttt{max.knot.dist}", "\\texttt{premade.knots}", "\\texttt{plot.knots}"),
    "Default Value" = c("\\texttt{'fourier'}", "\\texttt{8}", "\\texttt{'fourier'}", "\\texttt{6}", "\\texttt{4*diff(range(coords[,1]))}", "\\texttt{4*diff(range(coords[,2]))}", "\\texttt{floor(n.times/10)}", "\\texttt{n.times}", "\\texttt{2}", "\\texttt{mean(dist(coords))}", "\\texttt{NULL}", "\\texttt{FALSE}"),
    "Description" = c(
      "Indicates which type of basis functions to use for the linear and seasonal components. The default is \\texttt{'fourier'}. Other values accepted are \\texttt{'grid'} (for multiresolution bisquare bases) and \\texttt{'tps'} (for thin-plate splines).",
      "Number of basis functions to use for the linear and seasonal components. For Fourier bases, this value is $R_s$. For bisquare bases, this value is ignored. For thin-plate spline bases, the number of bases is \\texttt{floor(sqrt(n.spatial.bases))\\^2} to create an even grid.", 
      "The same as \\texttt{spatial.style} but for the factor loadings.  This does not have to be the same bases as \\texttt{spatial.style}.", 
      "The same as \\texttt{n.spatial.bases} but for the factor loadings.  This does not have to be the same vvalue as \\texttt{n.spatial.bases}.",
      "Frequency of the Fourier basis function for longitude (or, if using other coordinate system, the first coordinate value). This value is $f_{s_{[1]}}$. Default value is four times the range of the longitude coordinates. If not using Fourier bases, this argument is not used.", 
      "Same as \\texttt{freq.lon} but for latitude (or the second coordinate value).", 
      "Number of Fourier basis functions to use for the temporally-dependent factors. This value is $R_t$. The default value is \\texttt{floor(n.times/10)}.", 
      "Frequency of the Fourier bases functions for the temporally-dependent factors. This value is $f_t$. Default value is $T$ (\\texttt{n.times}).",
      "The number of resolutions when using the bisquare basis functions. If not using bisquare bases, this argument is not used.",
      "The distance beyond which a location is considered 'too far' from a knot, meaning its basis function value associated with that knot evaluates to zero. If not using bisquare bases, this argument is not used.",
      "A list of coordinates containing pre-specified knots. Each element of the list is a resolution. Each resolution should have the same number of columns as \\texttt{coords}. If not using bisquare bases, this argument is not used.",
      "TRUE/FALSE value indicating whether to provide a base R plot of the knot resolutions overlaid on top of the given \\texttt{coords}. If not using bisquare bases, this argument is not used."
    ),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Arguments to \\texttt{BSTFA} and \\texttt{BSTFAfull} associated with basis functions used for various components of the model.", 
  linesep = ""
) %>%  
  column_spec(3, width = "3in")

## ----dataplotsetup, echo=FALSE, warning=FALSE, include=FALSE, fig.width=3-----
map_data_loc <- ggplot2::map_data('state')[ggplot2::map_data('state')$region == 'utah',]
full_map <- ggplot2::map_data('state')
sf_polygon <- sf::st_sfc(sf::st_polygon(list(as.matrix(map_data_loc[,c(1,2)]))), crs=4326)
fixed_locations <- utahDataList$Coords[c(144,89,129,78),]

m = ggplot() +
  ## First layer: worldwide map
  geom_polygon(data = full_map,
               aes(x=long, y=lat, group = group),
               color = '#9c9c9c', fill = '#f3f3f3') +
  ## Second layer: Country map
  geom_polygon(data = map_data_loc,
               aes(x=long, y=lat, group = group),
               color = 'black', fill='#f3f3f3') +
  coord_map() +
  coord_fixed(1.3,
              xlim = c(min(utahDataList$Coords[,1])-0.5, max(utahDataList$Coords[,1])+0.5),
              ylim = c(min(utahDataList$Coords[,2])-0.5, max(utahDataList$Coords[,2])+0.5)) + 
  geom_point(data=utahDataList$Coords,
             aes(x=Longitude,y=Latitude),
             color='black', cex=1) +
  geom_point(data=fixed_locations,
             aes(x=Longitude,y=Latitude),
             color='red', cex=2.5) + 
  theme(axis.text=element_blank(),
        axis.ticks = element_blank()) + 
  xlab("") + 
  ylab("")

## ----fig-dataplot, echo=FALSE, fig.cap="Plot of Utah showing locations of fixed factors (in red) and all locations (in black).", fig.align="center", fig.width=3----
m

## ----fig-factorexample, fig.cap="Posterior mean (black line) and 95% credible interval (gray band) for the fourth factor.", fig.height=3.5, fig.width=6, fig.align="center"----
plot_factor(bstfa.output,
            together=FALSE,
            include.legend=FALSE,
            factor=4,
            type='mean')

## ----fig-mapexample, cache=FALSE, warning=FALSE, fig.cap="Posterior mean estimates across the state of Utah for the fourth loading.", fig.align="center", fig.width=3.5----
map_spatial_param(bstfa.output,
         parameter='loading',
         loading=4,
         type='mean',
         map=TRUE,
         state=TRUE,
         location='utah',
         fine=50)

## ----eval=FALSE---------------------------------------------------------------
# out <- BSTFA(ymat=utahDataList$TemperatureVals,
#       dates=utahDataList$Dates,
#       coords=utahDataList$Coords,
#       spatial.style='fourier',
#       load.style='fourier',
#       n.spatial.bases=8,
#       n.load.bases=6,
#       freq.lon=40,
#       freq.lat=30,
#       n.temp.bases=floor(nrow(utahDataList$TemperatureVals)/10),
#       freq.temp=nrow(utahDataList$TemperatureVals))

## ----fig-plottingknots, fig.cap="Location of the observed data (open black circles) and the bisquare bases knots for resolution 1 (red dots) and resolution 2 (green dots) for default knot locations.", fig.width=4, fig.height=4, fig.align="center"----
bstfa.plot_knots = BSTFA(ymat=utahDataList$TemperatureVals,
                         dates=utahDataList$Dates,
                         coords=utahDataList$Coords,
                         spatial.style='grid',
                         load.style='grid',
                         knot.levels=2,
                         plot.knots=TRUE,
                         verbose=FALSE,
                         iters=5)

## ----fig-customknots, fig.cap="Location of the observed data (open black circles) and the bisquare bases knots for resolution 1 (red dots) and resolution 2 (green dots) for user-defined knots.", fig.width=4, fig.height=4, fig.align="center"----
knots=list()
max.lon = max(utahDataList$Coords[,1])
min.lon = min(utahDataList$Coords[,1])
max.lat = max(utahDataList$Coords[,2])
min.lat = min(utahDataList$Coords[,2])
range.lon = max.lon-min.lon
range.lat = max.lat-min.lat
knots[[1]] = expand.grid(c(min.lon+(range.lon/4), min.lon+3*(range.lon/4)),
            c(min.lat+(range.lat/4), min.lat+3*(range.lat/4)))
knots[[2]] = expand.grid(c(min.lon+(range.lon/6), 
                           min.lon+(range.lon/2), 
                           min.lon+5*(range.lon/6)),
                         c(min.lat+(range.lat/6), 
                           min.lat+(range.lat/2), 
                           min.lat+5*(range.lat/6)))
bstfa.custom_knots = BSTFA(ymat=utahDataList$TemperatureVals,
                           dates=utahDataList$Dates,
                           coords=utahDataList$Coords,
                           spatial.style='grid',
                           load.style='grid',
                           knot.levels=2,
                           plot.knots=TRUE,
                           premade.knots=knots,
                           verbose=FALSE,
                           iters=5)

## ----eval=FALSE---------------------------------------------------------------
# bstfaTPS <- BSTFA(ymat=utahDataList$TemperatureVals,
#       dates=utahDataList$Dates,
#       coords=utahDataList$Coords,
#       spatial.style='tps',
#       load.style='tps',
#       n.spatial.bases=8,
#       n.load.bases=10)

## ----eval=FALSE---------------------------------------------------------------
# loglik = computeLogLik(bstfa.output,
#                        verbose=FALSE)
# loo::waic(loglik)
# loo::loo(loglik)

## ----tab-waic, results='asis', label="tab-waic", echo=FALSE-------------------
kable(
  data.frame(
    "Type" = rep(rep(c("Fourier", "Bisquare", "TPS"), each=2), 2),
    "\\# of Spatial Bases" = rep(c(8, 16, "4 (1-level)", "20 (2-levels)", 9, 16), 2),
    "\\# of Temporal Bases" = rep(c(126, 200), each=6), 
  "LOO-CV" = c(202.7, 203.2, 203.4, 203.6, "\\textbf{202.6}", 203.2, 204.6, 205.7, 204.7, 203.9, 205.2, 204.8),
  "WAIC" = c(195.8, 195.4, 195.4, 195.3, 195.7, 195.6, 191.6, "\\textbf{191.4}", 191.6, 191.9, 191.7, 191.8),
    check.names = FALSE
  ),
  format = 'latex',
  booktabs = TRUE,
  escape = FALSE,
  longtable = FALSE,
  caption = "Summary of model performance diagnostics, LOO-CV and WAIC, for different basis functions.",
  align="c",
  linesep = ""
)

## ----fig-basiscompare, fig.cap="Posterior mean estimates of the second loading for different style of spatial bases: Fourier (left), bisquare (center), and thin-plate spline (right).", out.width="6.5in", echo=F----
knitr::include_graphics("basis_comparison3.png")

## ----eval=F-------------------------------------------------------------------
# plot_trace(bstfa.output,
#            parameter='beta',
#            param.range=c(27),
#            density=FALSE)

## ----eval=FALSE---------------------------------------------------------------
# convergence_diag(bstfa.output,
#                   type='eSS',
#                   cutoff = 100)

