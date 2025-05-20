### Functions to perform posterior inference

### Plotting functions for STFA

#' Prediction
#' @param out Output from BSTFA or BSTFAfull.
#' @param location Either a single integer indicating the location in the data set to provide predictions or a vector of length 2 providing the longitude and latitude of the new location. If \code{location=NULL} (default), the function will return predictions for all in-sample locations.
#' @param type One of \code{all}, \code{mean} (default), \code{median}, \code{ub}, or \code{lb} indicating which summary statistic of the predicted values to return.
#' @param ci.level If \code{type='lb'} or \code{'ub'}, the percentiles for the posterior interval.
#' @param new_x If the original model included covariates \code{x}, include the same covariates for prediction \code{location}.
#' @param pred.int Logical scalar indicating whether to include additional uncertainty for posterior predictive intervals (\code{TRUE}; default) or not (posterior draws from the mean of \code{location}).
#' @returns A matrix or vector of predicted values for \code{location}.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' loc1means <- predictBSTFA(out, location=1, pred.int=FALSE)
#' @importFrom npreg basis.tps
#' @export predictBSTFA
predictBSTFA = function(out, location=NULL, type='mean',
                       ci.level = c(0.025, 0.975), new_x=NULL, pred.int=TRUE) {

  # FIX ME - do useful functions like bisquare still work?
  # FIX ME - add resid to "predict all locations"

  if (is.null(location)) { # predict for all observed locations
    facts <- matrix(0, ncol=out$draws, nrow=out$n.times*out$n.locs)
    for(i in 1:out$draws){
      facts[,i] <- c(matrix(out$F.tilde[i,],nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)%*%t(matrix(out$Lambda.tilde[i,],nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)))
    }
    ypreds = kronecker(Matrix::Diagonal(out$n.locs), rep(1,out$n.times))%*%t(out$mu) +
      kronecker(Matrix::Diagonal(out$n.locs), out$model.matrices$linear.Tsub)%*%t(out$beta) +
      kronecker(Matrix::Diagonal(out$n.locs), out$model.matrices$seasonal.bs.basis)%*%t(out$xi) +
      facts

  } else if (is.null(dim(location))) {  # predict a specific observed location
    loc.seq=c()
    xi.seq=c()
    lam.seq=c()
    for (i in 1:length(location)) {
      loc.seq <- append(loc.seq, ((location[i]-1)*out$n.times + 1):(location[i]*out$n.times))
      xi.seq <- append(xi.seq, ((location[i]-1)*out$n.seasn.knots + 1):(location[i]*out$n.seasn.knots))
      lam.seq <- append(lam.seq, ((location[i]-1)*out$n.factors + 1):(location[i]*out$n.factors))
    }

    facts <- matrix(0, ncol=out$draws, nrow=length(loc.seq))
    for(i in 1:out$draws){
      facts[,i] <- c(matrix(out$F.tilde[i,],nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)%*%t(matrix(out$Lambda.tilde[i,lam.seq],nrow=length(location),ncol=out$n.factors,byrow=TRUE)))
    }
    if (pred.int) {
      resid = matrix(rnorm(out$draws*out$n.times,
                           mean=rep(0,out$draws*out$n.times),
                           sd=sqrt(rep(c(out$sig2),each=out$n.times))),ncol=out$draws,byrow=TRUE)
      ypreds = kronecker(diag(length(location)), rep(1,out$n.times))%*%matrix(t(out$mu)[location,],nrow=length(location)) +
        kronecker(diag(length(location)), out$model.matrices$linear.Tsub)%*%matrix(t(out$beta)[location,],nrow=length(location)) +
        kronecker(diag(length(location)), out$model.matrices$seasonal.bs.basis)%*%t(out$xi)[xi.seq,] +
        facts + resid
    } else {
      ypreds = kronecker(diag(length(location)), rep(1,out$n.times))%*%matrix(t(out$mu)[location,],nrow=length(location)) +
        kronecker(diag(length(location)), out$model.matrices$linear.Tsub)%*%matrix(t(out$beta)[location,],nrow=length(location)) +
        kronecker(diag(length(location)), out$model.matrices$seasonal.bs.basis)%*%t(out$xi)[xi.seq,] +
        facts
    }

  } else if (length(dim(location))>1) { # predict at a new location (coordinates should have given to location)

    if (out$spatial.style=='grid') {
      # predS=makePredS(out,location)
      predS <- NULL
      for(kk in 1:length(out$knots.spatial)) {
        bspred <- bisquare2d(as.matrix(location), as.matrix(out$knots.spatial[[kk]]))
        predS <- cbind(predS, bspred)
      }
    }

    if (out$spatial.style == 'fourier') {
      m.fft.lon <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,1])/out$freq.lon)
        cos_term <- cos(2 * pi * k * (location[,1])/out$freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,2])/out$freq.lat)
        cos_term <- cos(2 * pi * k * (location[,2])/out$freq.lat)
        cbind(sin_term, cos_term)
      })
      Slon <- cbind(m.fft.lon[1:nrow(location),], m.fft.lon[(nrow(location)+1):(2*nrow(location)),])
      Slat <- cbind(m.fft.lat[1:nrow(location),], m.fft.lat[(nrow(location)+1):(2*nrow(location)),])
      predS = matrix(Slat*Slon,ncol=out$n.spatial.bases)
    }
    if (out$spatial.style == 'tps') {
      coords_added = rbind(out$coords,location)
      predS = matrix(npreg::basis.tps(coords_added, knots=out$knots.spatial, rk=TRUE)[-(1:nrow(out$coords)),-(1:2)],ncol=out$n.spatial.bases)
    }

    if (!is.null(new_x)) {
      predS <- cbind(predS, new_x)
      predloc <- predloc[complete.cases(predS),]
      predS <- predS[complete.cases(predS),]
    }

    ### Mu
    mumean <- predS%*%t(out$alpha.mu)
    if(pred.int){
        muresid = matrix(rnorm(nrow(location)*out$draws,
                             mean=rep(0,nrow(location)*out$draws),
                             sd=sqrt(rep(c(out$tau2.mu),each=nrow(location)))),ncol=out$draws,byrow=TRUE)
        mupred <- mumean + muresid
    }else{
     mupred <- mumean
    }
    mulong = kronecker(Matrix::Diagonal(nrow(location)),
                       rep(1,out$n.times))%*%mupred

    ### Beta (linear slope)
    betamean = predS%*%t(out$alpha.beta)
    if(pred.int){
        betaresid = matrix(rnorm(nrow(location)*out$draws,
                             mean=rep(0,nrow(location)*out$draws),
                             sd=sqrt(rep(c(out$tau2.beta),each=nrow(location)))),ncol=out$draws,byrow=TRUE)
        betapred <- betamean + betaresid
    }else{
        betapred <- betamean
    }
    betalong = kronecker(Matrix::Diagonal(nrow(location)),
                         out$model.matrices$linear.Tsub)%*%betapred

    ### Xi (seasonal)
    predS.xi = as(kronecker(predS, diag(out$n.seasn.knots)), "sparseMatrix")
    ximean <- predS.xi%*%t(out$alpha.xi)
    if(pred.int){
        xiresid <- matrix(rnorm(nrow(location)*out$n.seasn.knots*out$draws,
                            mean=rep(0,nrow(location)*out$n.seasn.knots*out$draws),
                            sd=sqrt(rep(c(out$tau2.xi),each=nrow(location)*out$n.seasn.knots))),ncol=out$draws,byrow=TRUE)
        xipred <- ximean + xiresid
    }else{
        xipred <- ximean
    }
    xilong = kronecker(Matrix::Diagonal(nrow(location)),
                       out$model.matrices$seasonal.bs.basis)%*%xipred

    ### Factor Analysis
    if (out$load.style == 'grid') {
      predQS <- NULL
      for(kk in 1:length(out$knots.load)) {
        bspred <- bisquare2d(as.matrix(location), as.matrix(out$knots.load[[kk]]))
        predQS <- cbind(predS, bspred)
      }
    }
    if (out$load.style == 'fourier') {
      m.fft.lon <- sapply(1:(out$n.load.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,1])/out$freq.lon)
        cos_term <- cos(2 * pi * k * (location[,1])/out$freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(out$n.load.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,2])/out$freq.lat)
        cos_term <- cos(2 * pi * k * (location[,2])/out$freq.lat)
        cbind(sin_term, cos_term)
      })
      Slon <- cbind(m.fft.lon[1:nrow(location),], m.fft.lon[(nrow(location)+1):(2*nrow(location)),])
      Slat <- cbind(m.fft.lat[1:nrow(location),], m.fft.lat[(nrow(location)+1):(2*nrow(location)),])
      predQS = matrix(Slat*Slon,ncol=out$n.load.bases)
    }
    if (out$load.style == 'tps') {
      coords_added = rbind(out$coords,location)
      predQS = matrix(npreg::basis.tps(coords_added, knots=out$knots.load, rk=TRUE)[-(1:nrow(out$coords)),-(1:2)],ncol=out$n.load.bases)
    }

    # Lambda (loadings)
    Lam = array(dim=c(nrow(predQS),out$n.factors,out$draws))
    for (i in 1:out$draws) {
      Lammean = predQS%*%matrix(out$alphaS[i,],nrow=out$n.load.bases,ncol=out$n.factors,byrow=TRUE)
      if(pred.int){
          Lamresid = matrix(rnorm(nrow(location)*out$n.factors,
                              mean=rep(0,nrow(location)*out$n.factors),
                              sd=sqrt(rep(c(out$tau2.lambda[i,]),each=out$n.factors))),ncol=out$n.factors,byrow=TRUE)
          Lam[,,i] = Lammean + Lamresid
      }else{
          Lam[,,i] = Lammean
      }
    }

    # F (factor scores)
    facts = array(dim=c(out$n.times,nrow(location),out$draws))
    for (i in 1:out$draws) {
      facts[,,i] = matrix(out$F.tilde[i,],nrow=out$n.times,ncol=out$n.factors)%*%matrix(t(Lam[,,i]),nrow=out$n.factors,ncol=nrow(location))
    }

    if (pred.int) {
      resid = matrix(rnorm(out$draws*out$n.times,
                           mean=rep(0,out$draws*out$n.times),
                           sd=sqrt(rep(c(out$sig2),each=out$n.times))),ncol=out$draws,byrow=TRUE)
      ypreds = mulong + betalong + xilong + matrix(facts, nrow=out$n.times*nrow(location), ncol=out$draws) + resid
    } else {
      ypreds = mulong + betalong + xilong + matrix(facts, nrow=out$n.times*nrow(location), ncol=out$draws)
    }

  }

  if (type == 'all') ypreds.return = ypreds

  if (type == 'mean') {
    ypreds.return = apply(ypreds,1,mean)
  }

  if (type == 'median') {
    ypreds.return = apply(ypreds,1,quantile,prob=c(0.5))
  }

  if (type == 'lb') {
    ypreds.return = apply(ypreds,1,quantile,prob=ci.level[1])
  }

  if (type == 'ub') {
    ypreds.return = apply(ypreds,1,quantile,prob=ci.level[2])
  }

  ypreds.return
}


#' Plot a location's time series of estimated/interpolated values.
#' @param out Output from BSTFA or BSTFAfull.
#' @param location Either a single integer indicating the location in the data set to plot or a vector of length 2 providing the longitude and latitude of the new location.
#' @param new_x If the original model included covariates \code{x}, include the same covariates for prediction \code{location}.
#' @param type One of \code{mean} (default), \code{median}, \code{ub}, or \code{lb} indicating which summary statistic of the predicted values to return.
#' @param par.mfrow A vector of length 2 indicating the number of rows and columns to divide the plotting window. Default is \code{c(1,1)}.
#' @param ci.level If \code{type='lb'} or \code{'ub'}, the percentiles for the posterior interval.
#' @param uncertainty Logical scalar indicating whether to plot the uncertainty bounds (\code{TRUE}; default) or not.
#' @param xrange A date vector of length 2 providing the lower and upper bounds of the dates to include in the plot.
#' @param truth Logical scalar indicating whether to include the observed measurements (\code{TRUE}) or not (default).  If \code{TRUE}, \code{location} must be an integer corresponding to the column of the data matrix for the in-sample prediction location.
#' @param ylim Numeric vector of length 2 providing the lower and upper bounds of the y-axis.  If \code{NULL} (default), the y-axis limits are chosen using the range of the predictions.
#' @returns A plot of predicted values for \code{location}.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.location(out, location=1, pred.int=FALSE)
#' @export plot.location
plot.location = function(out, location, new_x=NULL,
                         type='mean', par.mfrow=c(1,1), pred.int=TRUE,
                         ci.level = c(0.025, 0.975),
                         uncertainty=TRUE, xrange=NULL, truth=FALSE,
                         ylim=NULL) {

  ypreds = predictBSTFA(out, location=location, type=type, new_x=new_x, pred.int=pred.int)
  if (uncertainty) {
    ypreds.lb = predictBSTFA(out, location=location, type='lb',
                            new_x=new_x,
                            ci.level=ci.level,
                            pred.int=pred.int)
    ypreds.ub = predictBSTFA(out, location=location, type='ub',
                            new_x=new_x,
                            ci.level=ci.level,
                            pred.int=pred.int)
  }

  if (is.null(dim(location))) n.col=length(location)
  if (length(dim(location))>1) n.col=nrow(location)

  ymat.preds = matrix(ypreds, nrow=out$n.times, ncol=n.col)
  if (uncertainty) {
    ymat.preds.lb = matrix(ypreds.lb, nrow=out$n.times, ncol=n.col)
    ymat.preds.ub = matrix(ypreds.ub, nrow=out$n.times, ncol=n.col)
  }

  if (is.null(xrange)) xlims=1:out$n.times
  else xlims=which(out$dates > xrange[1] & out$dates < xrange[2])

  par(mfrow=par.mfrow)

  for (i in 1:n.col) {
    if (is.null(ylim)) {
      if (uncertainty) {
        if (truth) ylim = range(c(ymat.preds.ub[,i], ymat.preds.lb[,i], out$ymat[,location]),na.rm=TRUE)
        else ylim = range(c(ymat.preds.ub[,i], ymat.preds.lb[,i]))
      }
      else {
        if (truth) ylim = range(c(ymat.preds[,i], out$ymat[,location]),na.rm=TRUE)
        else ylim = range(ymat.preds[,i])
      }
    }
    plot(y=ymat.preds[xlims,i],
         x=out$dates[xlims],
         type='l',
         main = ifelse(is.null(dim(location)),
                       paste("Location", location[i]),
                       paste("Longitude", location[i,1], "Latitude", location[i,2])),
         xlab = "Time",
         ylab = "Value",
         ylim=ylim,
         lwd=2)
    mylegend <- c("Posterior Mean")
    mylines <- c(1)
    mydots <- c(NA)
    mylegcols <- c("black")
    mylwd <- c(1)
    if (uncertainty) {
      polygon(x=c(out$dates[xlims], rev(out$dates[xlims])), y=c(ymat.preds.lb[xlims,i], rev(ymat.preds.ub[xlims,i])), border=NA, col=rgb(.5, .5, .5, .4))
      #lines(ymat.preds.ub[xlims,i], x=out$dates[xlims], col='green', lty=2)
      mylegend <- c(mylegend, paste(100*(ci.level[2]-ci.level[1]), "% Uncertainty", sep=""))
      mylines <- c(mylines, 1)
      mydots <- c(mydots, NA)
      mylegcols <- c(mylegcols, rgb(.5, .5, .5, .4))
      mylwd <- c(mylwd, 3)
    }
    if (truth & is.null(dim(location))) {
      points(y=out$ymat[xlims,location[i]],x=out$dates[xlims], col=rgb(.5, .5, .5,1))
      mylegend <- c(mylegend, paste("Observations"))
      mylines <- c(mylines, NA)
      mydots <- c(mydots, 21)
      mylegcols <- c(mylegcols, rgb(.5, .5, .5, 1))
      mylwd <- c(mylwd, NA)
    }
    legend("topright", legend=mylegend, lty=mylines, lwd=mylwd, col=mylegcols, pch=mydots)
  }


}



#' Plot the spatially-dependent parameter for in-sample locations.
#' @param out Output from BSTFA or BSTFAfull.
#' @param parameter One of \code{"slope"} (default), \code{"loading"}, or \code{"mean"}.
#' @param loadings If \code{parameter="loading"}, an integer indicating which factor loading to plot.
#' @param type One of \code{mean} (default), \code{median}, \code{ub}, or \code{lb} indicating which summary statistic to plot at each location.
#' @param yearscsale If \code{parameter="slope"}, a logical scalar indicating whether to translate it to a yearly scale (\code{TRUE}; default).
#' @param ci.level If \code{type='lb'} or \code{'ub'}, the percentiles for the posterior interval.
#' @param color.gradient The color palette to use for the plot.  Default is \code{colorRampPalette(rev(RColorBrewer::brewer.pal(9, name='RdBu')))(50)}.
#' @returns A plot of spatially-dependent parameter values for the observed locations.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.grid(out, parameter="slope")
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export plot.grid
plot.grid = function(out, parameter, loadings=1, type='mean', ci.level=c(0.025, 0.975), yearscale=TRUE,
                     color.gradient=colorRampPalette(rev(RColorBrewer::brewer.pal(9, name='RdBu')))(50)) {

  if (parameter=='slope') {
    if (type=='mean') vals = apply(out$beta,2,mean)
    if (type=='median') vals = apply(out$beta,2,quantile,prob=0.5)
    if (type=='lb') vals = apply(out$beta,2,quantile,prob=ci.level[1])
    if (type=='ub') vals = apply(out$beta,2,quantile,prob=ci.level[2])
  } else if (parameter=='mean') {
    if (type=='mean') vals = apply(out$mu,2,mean)
    if (type=='median') vals = apply(out$mu,2,quantile,prob=0.5)
    if (type=='lb') vals = apply(out$mu,2,quantile,prob=ci.level[1])
    if (type=='ub') vals = apply(out$mu,2,quantile,prob=ci.level[2])
  } else if (parameter=='loading') {
    if (type=='mean') vals = matrix(apply(out$Lambda.tilde,2,mean),nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)
    if (type=='median') vals = matrix(apply(out$Lambda.tilde,2,median),nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)
    if (type=='lb') vals = matrix(apply(out$Lambda.tilde,2,quantile,prob=ci.level[1]),nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)
    if (type=='ub') vals = matrix(apply(out$Lambda.tilde,2,quantile,prob=ci.level[2]),nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)
  }

  if (parameter=='slope' & yearscale) {
    vals <- vals*365.25/(out$doy[2] - out$doy[1])
  }

  max_value = max(abs(min(vals)),abs(max(vals)))
  min_value = -max_value

  if (parameter == 'slope' | parameter == 'mean') {
    print(ggplot(mapping=aes(x=out$coords[,1], y=out$coords[,2],
                             color=vals)) +
            geom_point() +
            ggtitle(ifelse(parameter=='slope', "Slope", "Mean")) +
            xlab('Longitude') +
            ylab('Latitude') +
            labs(color = ifelse(parameter=='slope', "Slope", "Mean")) +
      scale_colour_gradientn(colors=color.gradient,
                             name='Slope', limits = c(min_value, max_value)))
  }
  if (parameter == 'loading') {
    for (i in loadings) {
      mm <- ggplot(mapping = aes(x = out$coords[,1], y = out$coords[,2])) +
        geom_point(aes(color = vals[,i])) +  # Main values
        geom_point(aes(x = out$coords[out$factors.fixed[i], 1],
                       y = out$coords[out$factors.fixed[i], 2],
                       shape = "Fixed Location"),
                   color = "black", size = 2, stroke = 1) +  # Fixed location as open dot
        xlab('Longitude') +
        ylab('Latitude') +
        ggtitle(paste("Loading ", i, ", Fixed Location ", out$factors.fixed[i], sep = "")) +
        scale_colour_gradientn(colors = color.gradient,
                               name = paste('Loading', loadings),
                               limits = c(min_value, max_value)) +
        scale_shape_manual(name = "", values = c("Fixed Location" = 1)) +  # shape 1 = open circle
        guides(color = guide_colorbar(order = 1),
               shape = guide_legend(order = 2))
      print(mm)
    }
  }
}


#' Plot a map of interpolated spatially-dependent parameter values.
#' @param out Output from BSTFA or BSTFAfull.
#' @param parameter One of \code{"slope"} (default), \code{"loading"}, or \code{"mean"}.
#' @param loadings If \code{parameter="loading"}, an integer indicating which factor loading to plot.
#' @param type One of \code{mean} (default), \code{median}, \code{ub}, or \code{lb} indicating which summary statistic to plot at each location.
#' @param yearscsale If \code{parameter="slope"}, a logical scalar indicating whether to translate it to a yearly scale (\code{TRUE}; default).
#' @param ci.level If \code{type='lb'} or \code{'ub'}, the percentiles for the posterior interval.
#' @param fine Integer specifying the number of grid points along both the longitude and latitude directions used to interpolate the parameter. The resulting interpolation grid will contain \code{fine*fine} total locations. If \code{map=TRUE}, \code{state=TRUE}, and \code{location} is specified, the grid will be clipped to the boundaries of the specified state, removing locations outside of it.
#' @param color.gradient The color palette to use for the plot.  Default is \code{colorRampPalette(rev(RColorBrewer::brewer.pal(9, name='RdBu')))(fine)}.
#' @param with.uncertainty Logical scalar indicating whether to include lower and upper credible interval bounds for the parameter.  Default is \code{FALSE}.
#' @param map Logical scalar indicating whether to include a map. Default value is \code{FALSE}.  If \code{TRUE}, \code{location} must be specified.
#' @param state Logical scalar used when \code{map=TRUE} indicating whether the \code{location} is a state in the United States (\code{TRUE}) or a country (\code{FALSE}).
#' @param location Name of region to include in the map.  Fed to \code{region} in the function \code{ggplot2::map_data}.
#' @param addthin Integer indicating the number of saved draws to thin.  Default is to not thin any \code{addthin=1}.  This can save time when the object is from \code{BSTFAfull} and \code{parameter='loading'}.
#' @returns A plot of spatially-dependent parameter values for a grid of interpolated locations.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.map(out, parameter="slope", map=T, state=T, location='utah', fine=50)
#' @importFrom npreg basis.tps
#' @importFrom sf st_sfc
#' @importFrom sf st_polygon
#' @importFrom sf st_point
#' @importFrom ggpubr ggarrange
#' @importFrom RColorBrewer brewer.pal
#' @import sf
#' @export plot.map
plot.map = function(out, parameter='slope', loadings=1, type='mean',
                    yearscale=TRUE, new_x=NULL,
                    ci.level=c(0.025, 0.975), fine=100,
                    color.gradient=colorRampPalette(rev(RColorBrewer::brewer.pal(9, name='RdBu')))(fine),
                    with.uncertainty=FALSE, map=FALSE, state=FALSE, location=NULL,
                    addthin=1) {

  # FIX ME - do functions like bisquare2d work in this function?

  if (map) {
    if (!state) {
      map_data_loc <- ggplot2::map_data('world')[ggplot2::map_data('world')$region == location,]
      full_map <- ggplot2::map_data('world')
    }
    if (state) {
      map_data_loc <- ggplot2::map_data('state')[ggplot2::map_data('state')$region == location,]
      full_map = ggplot2::map_data('state')
    }

    predloc <- expand.grid(seq(min(map_data_loc[,1]),
                               max(map_data_loc[,1]), length=fine),
                           seq(min(map_data_loc[,2]),
                               max(map_data_loc[,2]), length=fine))
  } else {
    predloc <- expand.grid(seq(min(out$coords[,1]),
                               max(out$coords[,1]), length=fine),
                           seq(min(out$coords[,2]),
                               max(out$coords[,2]), length=fine))
  }
  names(predloc) <- c("Lon", "Lat")

  if (parameter=='loading') {
    plot.title = paste('Loading', loadings)
    if (out$load.style=='grid') {
      predS <- NULL
      for(kk in 1:length(out$knots.load)) {
        bspred <- bisquare2d(as.matrix(predloc), as.matrix(out$knots.load[[kk]]))
        predS <- cbind(predS, bspred)
      }
    }
    if (out$load.style=='fourier') {
      ### Original Fourier Method
      m.fft.lon <- sapply(1:(out$n.load.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (predloc[,1])/out$freq.lon)
        cos_term <- cos(2 * pi * k * (predloc[,1])/out$freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(out$n.load.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (predloc[,2])/out$freq.lat)
        cos_term <- cos(2 * pi * k * (predloc[,2])/out$freq.lat)
        cbind(sin_term, cos_term)
      })
      Slon <- cbind(m.fft.lon[1:nrow(predloc),], m.fft.lon[(nrow(predloc)+1):(2*nrow(predloc)),])
      Slat <- cbind(m.fft.lat[1:nrow(predloc),], m.fft.lat[(nrow(predloc)+1):(2*nrow(predloc)),])
      predS = matrix(Slat*Slon,ncol=out$n.load.bases)
    }
    if (out$load.style=='tps') {
      predS = npreg::basis.tps(predloc,knots=out$knots.load,rk=TRUE)[,-(1:2)]
    }
    if(out$load.style=='full'){
      predS = diag(1, dim(predloc)[1])
    }
  }
  else {
    if (out$spatial.style=='grid') {
      plot.title = 'Slope'
      predS <- NULL
      for(kk in 1:length(out$knots.spatial)) {
        bspred <- bisquare2d(as.matrix(predloc), as.matrix(out$knots.spatial[[kk]]))
        predS <- cbind(predS, bspred)
      }
    }
    if (out$spatial.style=='fourier') {
      ### Original Fourier Method
      m.fft.lon <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (predloc[,1])/out$freq.lon)
        cos_term <- cos(2 * pi * k * (predloc[,1])/out$freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (predloc[,2])/out$freq.lat)
        cos_term <- cos(2 * pi * k * (predloc[,2])/out$freq.lat)
        cbind(sin_term, cos_term)
      })
      Slon <- cbind(m.fft.lon[1:nrow(predloc),], m.fft.lon[(nrow(predloc)+1):(2*nrow(predloc)),])
      Slat <- cbind(m.fft.lat[1:nrow(predloc),], m.fft.lat[(nrow(predloc)+1):(2*nrow(predloc)),])
      predS = matrix(Slat*Slon,ncol=out$n.spatial.bases)
   }
    if (out$spatial.style=='tps') {
      predS = npreg::basis.tps(predloc,knots=out$knots.spatial,rk=TRUE)[,-(1:2)]
    }
  }

  if (!is.null(new_x)) predS <- cbind(predS, new_x)
  predloc <- predloc[complete.cases(predS),]
  predS <- predS[complete.cases(predS),]

  if (parameter=='slope') {
    legend.name = 'Slope'
    betamean <- predS%*%t(out$alpha.beta[seq(1, out$draws, by=addthin),])
    betaresid <- matrix(rnorm(fine^2*floor(out$draws/addthin),
                             mean=rep(0,fine^2*floor(out$draws/addthin)),
                             sd=sqrt(rep(c(out$tau2.beta[seq(1, out$draws, by=addthin),]),each=fine^2))),ncol=floor(out$draws/addthin),byrow=TRUE)
    # betapred <- betamean + betaresid
    betapred <- betamean
    if (yearscale) {
      pred <- betapred*365.25/(out$doy[2] - out$doy[1])
    } else {
      pred <- betapred
    }
  }
  if (parameter=='mean') {
    legend.name = 'Mean'
    mumean <- predS%*%t(out$alpha.mu[seq(1, out$draws, by=addthin),])
    muresid <- matrix(rnorm(fine^2*floor(out$draws/addthin),
                              mean=rep(0,fine^2*floor(out$draws/addthin)),
                              sd=sqrt(rep(c(out$tau2.mu[seq(1, out$draws, by=addthin),]),each=fine^2))),ncol=floor(out$draws/addthin),byrow=TRUE)
    # pred <- mumean + muresid
    pred <- mumean
  }

  if (parameter=='loading') {
    legend.name = paste('Loading', loadings)
    if(out$load.style=="full"){
        names(out$coords) <- c("Lon", "Lat")
        npred <- dim(predloc)[1]
        predloc2 <- rbind(out$coords, predloc)
        preddist <- as.matrix(dist(predloc2))
        condinds <- 1:out$n.locs
        lammean <- matrix(0, nrow=floor(out$draws/addthin), ncol=npred)
        lamresid <- matrix(0, nrow=floor(out$draws/addthin), ncol=npred)
        mycount <- 0
        for(d in seq(1, out$draws, by=addthin)){
            mycount <- mycount + 1
            bigmat <- out$tau2.lambda[d,loadings]*exp(-preddist/out$phi.lambda[d,loadings])
            A <- bigmat[condinds, condinds]
            B <- bigmat[condinds, -condinds]
            C <- bigmat[-condinds, -condinds]

            L <- chol(A)
            LB <- forwardsolve(t(L), B)
            part1 <- t(backsolve(L, LB))

            condvar <- C - part1%*%B
            lammean[mycount,] <- part1%*%out$Lambda.tilde[d, seq(loadings, out$n.factors*out$n.locs, by=out$n.factors)] #((loading-1)*out$n.locs) + (1:out$n.locs)]
            cholC <- chol(condvar)
            lamresid[mycount,] <- as.numeric(cholC%*%rnorm(npred))
        }
        pred <- t(lammean)
    }else{
      lammean <- predS%*%t(out$alphaS)[seq(loadings,out$n.load.bases*out$n.factors,by=out$n.factors),seq(1, out$draws, by=addthin)]
      lamresid <- matrix(rnorm(fine^2*floor(out$draws/addthin),
                             mean=rep(0,fine^2*floor(out$draws/addthin)),
                             sd=sqrt(rep(c(out$tau2.lambda[seq(1, out$draws, by=addthin),loadings]),each=fine^2))),ncol=floor(out$draws/addthin),byrow=TRUE)
      pred <- lammean
    }
  }

  if (type=='mean') {
    predloc$predm <- apply(pred, 1, mean)
    plot.title = paste0(toupper(substring(type,1,1)), substring(type,2))
  }
  if (type=='median') {
    predloc$predm <- apply(pred, 1, median)
    plot.title = paste0(toupper(substring(type,1,1)), substring(type,2))
  }
  if (type=='lb') {
    predloc$predm <- apply(pred, 1, quantile, prob=ci.level[1])
    plot.title = paste0((ci.level[2] - ci.level[1])*100, "% Lower Bound")
  }
  if (type=='ub') {
    predloc$predm <- apply(pred, 1, quantile, prob=ci.level[2])
    plot.title = paste0((ci.level[2] - ci.level[1])*100, "% Upper Bound")
  }
  if (!with.uncertainty) {
    max_value = max(abs(min(predloc$predm)),abs(max(predloc$predm)))
    min_value = -max_value
  }
  if (with.uncertainty) {
    predloc$predl <- apply(pred, 1, quantile, prob=ci.level[1])
    predloc$predu <- apply(pred, 1, quantile, prob=ci.level[2])
    max_value = max(abs(min(predloc$predl)),abs(max(predloc$predl)), abs(min(predloc$predu)),abs(max(predloc$predu)))
    min_value = -max_value
  }

  if (!map) {
    m <- ggplot2::ggplot(aes(x=Lon, y=Lat, fill=predm), data=predloc) +
      #geom_point(aes(x=Lon, y=Lat, color=predm)) +
      geom_raster(interpolate=T) +
      scale_fill_gradientn(colours=color.gradient, name=legend.name,
                             limits = c(min_value, max_value)) +
      ggtitle(plot.title) + xlab("Longitude") + ylab("Latitude")
    if (!with.uncertainty) print(m)
    if (with.uncertainty) {
      l <- ggplot(data=predloc, aes(x=Lon, y=Lat, fill=predl)) +
        #geom_point(aes(x=Lon, y=Lat, color=predl)) +
        geom_raster(interpolate=T) +
        scale_fill_gradientn(colours=color.gradient, name=legend.name,
                               limits = c(min_value, max_value)) +
        ggtitle(paste0((ci.level[2]-ci.level[1])*100,'% Lower Bound')) + xlab("Longitude") + ylab("Latitude")
      u <- ggplot(data=predloc, aes(x=Lon, y=Lat, fill=predu)) +
        #geom_point(aes(x=Lon, y=Lat, color=predu)) +
        geom_raster(interpolate=T) +
        scale_fill_gradientn(colours=color.gradient, name=legend.name,
                               limits = c(min_value, max_value)) +
        ggtitle(paste0((ci.level[2]-ci.level[1])*100,'% Upper Bound')) + xlab("Longitude") + ylab("Latitude")
      print(ggpubr::ggarrange(l, m, u, nrow=1, common.legend=T, legend="right"))
    }
  }

  if (map) {

    sf_polygon <- sf::st_sfc(sf::st_polygon(list(as.matrix(map_data_loc[,c(1,2)]))), crs=4326)

    sf_polygon <- sf::st_make_valid(sf_polygon)
    if(!sf::st_is_valid(sf_polygon)){suppressMessages(sf::sf_use_s2(FALSE))}
    ### Check if points fall inside of polygon ###
    inside = c()
    for (kk in 1:nrow(predloc)) {
      point = sf::st_sfc(sf::st_point(as.matrix(predloc[kk,c(1,2)])), crs=4326)
      if (suppressMessages(sf::st_intersects(point, sf_polygon, sparse=FALSE))) inside = append(inside, kk)
    }

    predloc.inside = predloc[inside, ]
    if (!with.uncertainty) {
      max_value = max(abs(min(predloc.inside$predm)),abs(max(predloc.inside$predm)))
      min_value = -max_value
    } else{
      max_value = max(abs(min(predloc.inside$predl)),abs(max(predloc.inside$predl)), abs(min(predloc.inside$predu)),abs(max(predloc.inside$predu)))
      min_value = -max_value
    }

    m = ggplot() +
      ## First layer: worldwide map
      geom_polygon(data = full_map,
                   aes(x=long, y=lat, group = group),
                   color = '#9c9c9c', fill = '#f3f3f3') +
      ## Second layer: Country map
      geom_polygon(data = map_data_loc,
                   aes(x=long, y=lat, group = group),
                   color = '#9c9c9c', fill='#f3f3f3') +
      coord_map() +
      coord_fixed(1.3,
                  xlim = c(min(out$coords[,1])-1, max(out$coords[,1])+1),
                  ylim = c(min(out$coords[,2])-1, max(out$coords[,2])+1)) +
      ggtitle(plot.title) + # FIX ME
      theme(panel.background =element_rect(fill = rgb(0.67, .84, .89, .35))) +
      #geom_point(data=predloc.inside, aes(x=Lon, y=Lat, color=predm)) +
      geom_raster(data=predloc.inside, aes(x=Lon, y=Lat, fill=predm)) +
      scale_fill_gradientn(colours=color.gradient, name=legend.name,
                             limits = c(min_value, max_value)) +
                             # limits = c(-0.7, 0.7)) + # FIX ME
      xlab('Longitude') +
      ylab('Latitude')
    if(!with.uncertainty){print(m)}


    if (with.uncertainty) {
      l = ggplot() +
        ## First layer: worldwide map
        geom_polygon(data = full_map,
                     aes(x=long, y=lat, group = group),
                     color = '#9c9c9c', fill = '#f3f3f3') +
        ## Second layer: Country map
        geom_polygon(data = map_data_loc,
                     aes(x=long, y=lat, group = group),
                     color = 'gray', fill='gray') +
        coord_map() +
        coord_fixed(1.3,
                    xlim = c(min(out$coords[,1])-1, max(out$coords[,1])+1),
                    ylim = c(min(out$coords[,2])-1, max(out$coords[,2])+1)) +
        ggtitle(paste0((ci.level[2]-ci.level[1])*100,'% Lower Bound')) +
        theme(panel.background =element_rect(fill = rgb(0.67, .84, .89, .35))) +
        #geom_point(data=predloc.inside, aes(x=Lon, y=Lat, color=predl)) +
        geom_raster(data=predloc.inside, aes(x=Lon, y=Lat, fill=predl), interpolate=T) +
        scale_fill_gradientn(colours=color.gradient, name=legend.name,
                               limits = c(min_value, max_value)) +
        xlab('Longitude') +
        ylab('Latitude')

      u = ggplot() +
        ## First layer: worldwide map
        geom_polygon(data = full_map,
                     aes(x=long, y=lat, group = group),
                     color = '#9c9c9c', fill = '#f3f3f3') +
        ## Second layer: Country map
        geom_polygon(data = map_data_loc,
                     aes(x=long, y=lat, group = group),
                     color = 'gray', fill='gray') +
        coord_map() +
        coord_fixed(1.3,
                    xlim = c(min(out$coords[,1])-1, max(out$coords[,1])+1),
                    ylim = c(min(out$coords[,2])-1, max(out$coords[,2])+1)) +
        ggtitle(paste0((ci.level[2]-ci.level[1])*100,'% Upper Bound')) +
        theme(panel.background =element_rect(fill = rgb(0.67, .84, .89, .35))) +
        #geom_point(data=predloc.inside, aes(x=Lon, y=Lat, color=predu)) +
        geom_raster(data=predloc.inside, aes(x=Lon, y=Lat, fill=predu), interpolate=T) +
        scale_fill_gradientn(colours=color.gradient, name=legend.name,
                               limits = c(min_value, max_value)) +
        xlab('Longitude') +
        ylab('Latitude')

      print(ggpubr::ggarrange(l, m, u, nrow=1, common.legend=T, legend="right"))
    }
  }

}


#' Plot the factors
#' @param out Output from BSTFA or BSTFAfull.
#' @param factor Integer or vector of integers specifying which factor(s) to plot.
#' @param together If \code{length(factor)>1}, logical scalar specifying whether to plot all factors on a single plot. Default is \code{FALSE}.
#' @param include.legend If \code{length(factor)>1} and \code{together=TRUE}, a logical scalar specifying whether to include a legend.  Default is \code{TRUE}.
#' @param type One of \code{mean} (default), \code{median}, \code{ub}, or \code{lb} indicating which summary statistic to plot at each location.
#' @param uncertainty Logical scalar indicating whether to include lower and upper credible interval bounds for the parameter.  Default is \code{FALSE}.
#' @param ci.level A vector of length 2 specifying the quantiles to use for lower and upper bounds for \code{type='lb'}, \code{type='ub'}, or \code{uncertainty=TRUE}.
#' @param xrange A date vector of length 2 providing the lower and upper bounds of the dates to include in the plot.
#' @returns A plot of spatially-dependent parameter values for a grid of interpolated locations.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.factor(out, factor=1:4, together=T)
#' @export plot.factor
plot.factor = function(out, factor=1, together=FALSE, include.legend=TRUE,
                       type='mean', uncertainty=T, ci.level=c(0.025, 0.975),
                       xrange=NULL) {

  par(mfrow=c(1,1))

  if (type=='mean') F.tilde = matrix(apply(out$F.tilde,2,mean),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
  if (type=='median') F.tilde = matrix(apply(out$F.tilde,2,median),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
  if (type=='lb') F.tilde = matrix(apply(out$F.tilde,2,quantile,prob=ci.level[1]),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
  if (type=='ub') F.tilde = matrix(apply(out$F.tilde,2,quantile,prob=ci.level[2]),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
  if(uncertainty){
    F.tilde.lb <- matrix(apply(out$F.tilde,2,quantile,prob=ci.level[1]),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
    F.tilde.ub <- matrix(apply(out$F.tilde,2,quantile,prob=ci.level[2]),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
  }

  if (is.null(xrange)) xlims=1:out$n.times
  else xlims=which(out$dates > xrange[1] & out$dates < xrange[2])
  if (together) {
    mycols <- RColorBrewer::brewer.pal(out$n.factors, 'Dark2')
    mycolssee <- paste0(mycols, "20")
    plot(y=F.tilde[xlims,1], x=out$dates[xlims], type='l', main = ('All Factors'),
         xlab = 'Time', ylab='Value', col=mycols[1],
         ylim=range(F.tilde))
         # ylim=c(-7,7)) # FIX ME
    if(uncertainty){
      for(i in 1:out$n.factors){
      polygon(x=c(out$dates[xlims], rev(out$dates[xlims])), y=c(F.tilde.lb[xlims,i], rev(F.tilde.ub[xlims,i])), col=mycolssee[i], border=NA)
      }
    }
    for (i in 1:out$n.factors) {
      lines(y=F.tilde[,i], x=out$dates, type='l', col=mycols[i])
    }
    if (include.legend) {
      legend("topleft",
             legend=paste("Factor", seq(1,out$n.factors)),
             col = mycols,
             lty=1,
             lwd=2,
             xpd=TRUE)
             #inset=c(0.2,0))
    }

  }
  if (!together) {
    for (i in factor) {
      if(uncertainty){
        ylims <- range(c(c(F.tilde.lb, F.tilde.ub)))
      }else{ylims <- range(F.tilde)}
      plot(y=F.tilde[xlims,i], x=out$dates[xlims], type='l', main=paste('Factor', i),
           xlab = 'Time', ylab='Value', lwd=2, ylim=ylims)
      if(uncertainty){polygon(x=c(out$dates[xlims], rev(out$dates[xlims])), y=c(F.tilde.lb[xlims,i], rev(F.tilde.ub[xlims,i])), col=rgb(.5, .5, .5,.4), border=NA)}
    }
  }

}


#' Plot annual curve
#' @param out Output from BSTFA or BSTFAfull.
#' @param location Either a single integer indicating the location in the data set to plot or a vector of length 2 providing the longitude and latitude of the new location.
#' @param add Logical scalar indicating whether the annual/seasonal process should be added to the existing plot.  Default is \code{FALSE}.
#' @param years Either \code{"one"} (indicating to plot just a single year; default) or \code{"all"} (indicating to plot all years in the observed time period).
#' @param new_x If the original model included covariates \code{x}, include the same covariates for \code{location}.
#' @param interval Numeric value between 0 and 1 specifying the probability of the credible interval.
#' @param yrange Numeric vector of length 2 providing the lower and upper bounds of the y-axis.  If \code{NULL} (default), the y-axis limits are chosen using the range of the seasonal process and data.
#' @returns A plot of the annual/seasonal process at \code{location}.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.annual(out, location=1)
#' @importFrom mgcv cSplineDes
#' @export plot.annual
plot.annual <- function(out, location, add=F,
                        years="one",
                        interval=0.95, yrange=NULL,
                        new_x=NULL){

  y = out$y
  x_set = out$doy
  if(years=="all"){
    dates.pred <- seq(as.Date(out$dates[1]), as.Date(out$dates[length(out$dates)]), by=1)
    doy.pred <- as.numeric(strftime(dates.pred, format="%j"))
  }else{
    dates.pred <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by=1) # 2001 doesn't matter; just gets correct doy
    doy.pred <- as.numeric(strftime(dates.pred, format="%j"))
    months.plot <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by="month") # 2001 doesn't matter; just gets correct month
    at.doy.plot <- as.numeric(strftime(months.plot, format="%j"))
    months.plot <- months(months.plot, abbreviate=T)
  }

  knots <- seq(1, 366, length=out$n.seasn.knots+1)
  bs.basis <- mgcv::cSplineDes(doy.pred, knots)

  if(length(location)==1){
    loc.seq <- ((location-1)*out$n.seasn.knots + 1):(location*out$n.seasn.knots)
    xi.pred <- out$xi[,loc.seq]
    ann.pred <- bs.basis%*%t(xi.pred)
    ann.pred.mean <- apply(ann.pred, 1, mean)
    if(interval>0){
      ann.pred.bounds <- apply(ann.pred, 1, quantile, probs=c((1-interval)/2, (1+interval)/2))
    }else{
      ann.pred.bounds <- NULL
    }
    if(add==T){
      lines(dates.pred, ann.pred.mean, lwd=1.5)
      if(interval>0){
        polygon(c(dates.pred, rev(dates.pred)), c(ann.pred.bounds[1,], rev(ann.pred.bounds[2,])), col=rgb(.5, .5, .5, .4), border=NA)
      }
      lines(dates.pred, ann.pred.mean, lwd=2)
    }else{ #end if add==T
      if(length(y)>0){
        y.this <- out$y[((location-1)*out$n.times +1):(location*out$n.times)]
      }else{
        y.this <- NULL
      }
      if(is.null(yrange)==T){
        ylims <- range(c(ann.pred.mean, ann.pred.bounds, y.this), na.rm=T)
      }else{
        ylims <- yrange
      }
      if(years=="all"){
        plot(dates.pred, ann.pred.mean, lwd=1.5, type='l', xlab="Date", ylab="Annual Seasonal Cycle", ylim=ylims, main=paste("Location", location))
      }else{
        plot(doy.pred, ann.pred.mean, lwd=1.5, type='l', xaxt="n", xlab="Date", ylab="Annual Seasonal Cycle", ylim=ylims)
        axis(1, at=at.doy.plot, labels=months.plot)
      }
      if(interval>0){
        if(years=="all"){
          polygon(c(dates.pred, rev(dates.pred)), c(ann.pred.bounds[1,], rev(ann.pred.bounds[2,])), col=rgb(.5, .5, .5, .4), border=NA)
          lines(dates.pred, ann.pred.mean, lwd=1.5)
        }else{
          polygon(c(doy.pred, rev(doy.pred)), c(ann.pred.bounds[1,], rev(ann.pred.bounds[2,])), col=rgb(.5, .5, .5, .4), border=NA)
          lines(doy.pred, ann.pred.mean, lwd=1.5)
        }
      }

      if(length(y)>0){
        if(years=="all"){
          dates.data <- as.Date(out$dates)
          points(dates.data, y.this, col=rgb(.5, .5, .5, .25))
        }else{
          points(x_set, y.this, col=rgb(.5, .5, .5,.25))
        }
      }
    }
  } else { # New location

    if (out$spatial.style=='grid') {
      # predS=makePredS(out,location)
      predS <- NULL
      for(kk in 1:length(out$knots.spatial)) {
        bspred <- bisquare2d(as.matrix(location), as.matrix(out$knots.spatial[[kk]]))
        predS <- cbind(predS, bspred)
      }
    }

    if (out$spatial.style == 'fourier') {
      m.fft.lon <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,1])/out$freq.lon)
        cos_term <- cos(2 * pi * k * (location[,1])/out$freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(out$n.spatial.bases/2), function(k) {
        sin_term <- sin(2 * pi * k * (location[,2])/out$freq.lat)
        cos_term <- cos(2 * pi * k * (location[,2])/out$freq.lat)
        cbind(sin_term, cos_term)
      })
      Slon <- cbind(m.fft.lon[1:nrow(location),], m.fft.lon[(nrow(location)+1):(2*nrow(location)),])
      Slat <- cbind(m.fft.lat[1:nrow(location),], m.fft.lat[(nrow(location)+1):(2*nrow(location)),])
      predS = matrix(Slat*Slon,ncol=out$n.spatial.bases)
    }
    if (out$spatial.style == 'tps') {
      coords_added = rbind(out$coords,location)
      predS = matrix(npreg::basis.tps(coords_added, knots=out$knots.spatial, rk=TRUE)[-(1:nrow(out$coords)),-(1:2)],ncol=out$n.spatial.bases)
    }

    if (!is.null(new_x)) {
      predS <- cbind(predS, new_x)
      predloc <- predloc[complete.cases(predS),]
      predS <- predS[complete.cases(predS),]
    }

    predS.xi = as(kronecker(predS, diag(out$n.seasn.knots)), "sparseMatrix")
    ximean <- predS.xi%*%t(out$alpha.xi)
    xiresid <- matrix(rnorm(nrow(location)*out$n.seasn.knots*out$draws,
                            mean=rep(0,nrow(location)*out$n.seasn.knots*out$draws),
                            sd=sqrt(rep(c(out$tau2.xi),each=nrow(location)*out$n.seasn.knots))),ncol=out$draws,byrow=TRUE)
    # xi.pred <- ximean + xiresid
    xi.pred <- ximean
    ann.pred <- bs.basis%*%xi.pred
    ann.pred.mean <- apply(ann.pred, 1, mean)
    if(interval>0){
      ann.pred.bounds <- apply(ann.pred, 1, quantile, probs=c((1-interval)/2, (1+interval)/2))
    }else{
      ann.pred.bounds <- NULL
    }

    if(is.null(yrange)==T){
      ylims <- range(c(ann.pred.mean, ann.pred.bounds), na.rm=T)
    }else{
      ylims <- yrange
    }
    if(years=="all"){
      plot(dates.pred, ann.pred.mean, lwd=1.5, type='l', xlab="Date", ylab="Annual Seasonal Cycle", ylim=ylims, main=paste("Location", location[1], location[2]))
    }else{
      plot(doy.pred, ann.pred.mean, lwd=1.5, type='l', xaxt="n", xlab="Date", ylab="Annual Seasonal Cycle", ylim=ylims)
      axis(1, at=at.doy.plot, labels=months.plot)
    }
    if(interval>0){
      if(years=="all"){
        polygon(c(dates.pred, rev(dates.pred)), c(ann.pred.bounds[1,], rev(ann.pred.bounds[2,])), col=rgb(.5, .5, .5, .4), border=NA)
        lines(dates.pred, ann.pred.mean, lwd=1.5)
      }else{
        polygon(c(doy.pred, rev(doy.pred)), c(ann.pred.bounds[1,], rev(ann.pred.bounds[2,])), col=rgb(.5, .5, .5, .4), border=NA)
        lines(doy.pred, ann.pred.mean, lwd=1.5)
      }
    }
  }
}
