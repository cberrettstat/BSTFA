## ----setup1, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache=FALSE
)

## ----setup2, echo=FALSE, include=FALSE----------------------------------------
library(knitr)
devtools::load_all()
# library(BSTFA)

## ----fig-utahTemps, echo=FALSE, fig.height=4, fig.width=10, fig.cap='Mean-centered 30-day average daily minimum temperatures for 3 Utah weather stations (Moab, left; Canyonlands, middle; Logan, right) from 1999 to 2001.'----

colnames(utahDataList$TemperatureVals) = utahDataList$Locations
par(mfrow=c(1,3))
window=1050:1090
plot(y=utahDataList$TemperatureVals[window,which(utahDataList$Locations=='MOAB')],
     x=utahDataList$Dates[window], type = 'l',
     main = "Moab",
     xlab = "",
     ylab = "",
     ylim=c(-18,23),
     cex.main=1.5)
plot(y=utahDataList$TemperatureVals[window,which(utahDataList$Locations=='CANYONLANDS.THE.NECK')],
     x=utahDataList$Dates[window], type = 'l',
     main = "Canyonlands",
     xlab = "",
     ylab = "",
     ylim=c(-18,23),
     cex.main=1.5)
plot(y=utahDataList$TemperatureVals[window,utahDataList$Locations=='LOGAN.UTAH.ST.UNIV'], 
     x=utahDataList$Dates[window], type = 'l',
     main = "Logan",
     xlab = "",
     ylab = "",
     ylim=c(-18,23),
     cex.main=1.5)

## ----BSTFAFunction, echo=FALSE, include=FALSE, cache=FALSE--------------------
# This code runs the function correctly but isn't shown in the vignette
set.seed(23416)
bstfa.output = BSTFA(ymat=utahDataList$TemperatureVals,
                     dates=utahDataList$Dates,
                     coords=utahDataList$Coords,
                     factors.fixed = c(144,89,129,78), # Specific fixed factor locations for Utah data set
                     iters=1000,
                     burn=10,
                     thin=2,
                     verbose=FALSE,
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

## ----plotpredall--------------------------------------------------------------
loc = 1 # Alpine, Utah in our data set
preds = predictBSTFA(bstfa.output, 
                     location = loc,
                     type='all',
                     pred.int=TRUE,
                     ci.level=c(0.025,0.975))
dim(preds)

## ----plotpredobs--------------------------------------------------------------
loc = 1 # Alpine, Utah in our data set
preds = predictBSTFA(bstfa.output, 
                     location = loc,
                     type='mean',
                     pred.int=TRUE,
                     ci.level=c(0.025,0.975))
dim(preds)
length(preds)

## ----plotpred-----------------------------------------------------------------
loc = data.frame('Longitude' = 111.41, 'Latitude' = 38.29) # Torrey, Utah
preds_new = predictBSTFA(bstfa.output, 
                         location = loc,
                         type='all',
                         pred.int=TRUE,
                         ci.level=c(0.025,0.975))

## ----plotobserved-------------------------------------------------------------
loc = 1 # Alpine, Utah in our data set
plot_location(bstfa.output,
              location=loc,
              type='mean',
              pred.int=TRUE,
              uncertainty=TRUE,
              ci.level=c(0.025,0.975),
              xrange=c('1959-01-01', '1979-01-01'),
              truth=TRUE)

## ----plotprediction-----------------------------------------------------------
loc = data.frame('Longitude' = 111.41, 'Latitude' = 38.29) # Torrey, Utah
plot_location(bstfa.output,
              location=loc,
              type='mean',
              pred.int=TRUE,
              uncertainty=TRUE,
              ci.level=c(0.025,0.975),
              xrange=c('1959-01-01', '1979-01-01'),
              truth=FALSE)

## ----plotseason---------------------------------------------------------------
plot_annual(bstfa.output,
            location=1,
            years='one')

## ----plotobsslope-------------------------------------------------------------
plot_spatial_param(bstfa.output,
          type='mean',
          parameter='slope',
          yearscale=TRUE)

## ----plotobsloading1----------------------------------------------------------
plot_spatial_param(bstfa.output,
          type='mean',
          parameter='loading',
          loadings=1)

## ----factorexample2-----------------------------------------------------------
plot_factor(bstfa.output,
            together=FALSE,
            include.legend=FALSE,
            factor=1,
            type='mean')

## ----factorexample1-----------------------------------------------------------
plot_factor(bstfa.output,
            together=TRUE,
            include.legend=TRUE,
            type='mean')

## ----mapexample1, cache=FALSE, warning=FALSE----------------------------------
map_spatial_param(bstfa.output,
         parameter='slope',
         yearscale=TRUE,
         type='mean',
         map=TRUE,
         state=TRUE,
         location='utah',
         fine=100)

## ----comptime-----------------------------------------------------------------
compute_summary(bstfa.output)

## ----fourierplots-------------------------------------------------------------
plot_fourier_bases(utahDataList$Coords,
                   R=6,
                   plot.3d=TRUE,
                   freq.lon = diff(range(utahDataList$Coords[,1]))^2,
                   freq.lat = diff(range(utahDataList$Coords[,2]))^2)

## ----dataplotsetup, echo=FALSE, warning=FALSE, include=FALSE------------------
library(ggplot2)
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
             color='black', cex=0.5) +
  geom_point(data=fixed_locations,
             aes(x=Longitude,y=Latitude),
             color='red', cex=5) + 
  theme(axis.text=element_blank(),
        axis.ticks = element_blank()) + 
  xlab("") + 
  ylab("")

## ----dataplot, echo=FALSE-----------------------------------------------------
m

## ----factorexample------------------------------------------------------------
plot_factor(bstfa.output,
            together=FALSE,
            include.legend=FALSE,
            factor=4,
            type='mean')

## ----mapexample, cache=FALSE, warning=FALSE-----------------------------------
map_spatial_param(bstfa.output,
         parameter='loading',
         loading=4,
         yearscale=TRUE,
         type='mean',
         map=TRUE,
         state=TRUE,
         location='utah',
         fine=100)

## ----eval=FALSE---------------------------------------------------------------
# BSTFA(ymat=utahDataList$TemperatureVals,
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

## ----plottingknots------------------------------------------------------------
bstfa.plot_knots = BSTFA(ymat=utahDataList$TemperatureVals,
                         dates=utahDataList$Dates,
                         coords=utahDataList$Coords,
                         spatial.style='grid',
                         load.style='grid',
                         knot.levels=2,
                         plot.knots=TRUE,
                         verbose=FALSE,
                         iters=10)

## ----custom.knots-------------------------------------------------------------
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
                           iters=10)

## ----eval=FALSE---------------------------------------------------------------
# BSTFA(ymat=utahDataList$TemperatureVals,
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

## -----------------------------------------------------------------------------
plot_trace(bstfa.output,
           parameter='beta',
           param.range=c(27),
           density=FALSE)

## ----eval=FALSE---------------------------------------------------------------
# convergence_diag(bstfa.output,
#                   type='eSS',
#                   cutoff = 100)

