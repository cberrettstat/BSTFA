##### Functions for bstfa objects #####

#' Summary for 'bstfa' Objects
#'
#' This function provides posterior mean, standard deviations, and probabilities of being greater than 0 of interpretable parameters. 
#' @param object A bstfa object.
#' @param ... Additional arguments passed to or from other methods.
#' @returns A data frame containing the following elements:
#' \describe{
#'   \item{post.mean}{Posterior mean of parameters for each observed location.}
#'   \item{post.sd}{Posterior standard deviation of parameters for each observed location.}
#'   \item{post.pgt0}{Posterior probability that the parameter is greater than 0, if relevant.}
#' }
#' @author Candace Berrett and Adam Simpson
#' @examples
#' summary(out.sim)
#' @export 
summary.bstfa <- function(object, ...){
  
  #Collect row names for data frame
  station.names <- row.names(object$coords)
  
  #posterior summaries of spatially-dependent mean
  if(object$mean==TRUE){
    postm.mu <- apply(object$mu, 2, mean)
    postsd.mu <- apply(object$mu, 2, sd)
    postp.mu <- apply(object$mu>0, 2, mean)
    df.names <- paste0(station.names, ".mu")
  }else{
    postm.mu <- NULL
    postsd.mu <- NULL
    postp.mu <- NULL
    df.names <- NULL
  }
  
  #posterior summaries of spatially-dependent linear component
  if(object$linear==TRUE){
    postm.beta <- apply(object$beta, 2, mean)
    postsd.beta <- apply(object$beta, 2, sd)
    postp.beta <- apply(object$beta>0, 2, mean)
    df.names <- c(df.names, paste0(station.names, ".beta"))
  }else{
    postm.beta <- NULL
    postsd.beta <- NULL
    postp.beta <- NULL
    df.names <- c(df.names, NULL)
  }
  
  #posterior summary of residual standard deviation
  postm.sig <- mean(sqrt(object$sig2))
  postsd.sig <- sd(sqrt(object$sig2))
  postp.sig <- NA
  df.names <- c(df.names, 'sigma')
  
  mysum <- data.frame(post.mean=c(postm.mu, postm.beta, postm.sig), post.sd=c(postsd.mu, postsd.beta, postsd.sig), post.pgt0=c(postp.mu, postp.beta, postp.sig))
  row.names(mysum) <- df.names
  
  return(mysum)
  
}


#' Basic Plotting for 'bstfa' Objects
#'
#' This function calls appropriate plotting functions for different bstfa parameters. For more control over specific plotting features, use \code{plot_spatial_param}, \code{plot_factor}, \code{plot_annual}, \code{plot_location}, and \code{map_spatial_param}. 
#' @param x A bstfa object.
#' @param ... Additional arguments passed to or from other methods.
#' @param plot.param Parameter to plot: one of 'mu', 'beta', 'lambda', 'xi', 'annual', or 'factors'.
#' @param type One of 'map' (default; see \code{map_spatial_param}), 'basic' (see \code{plot_spatial_param}), 'location' (see \code{plot_location}), 'annual' (see \code{plot_annual}), or 'factor' (see \code{plot_factor}).
#' @param location If \code{type='annual'} or \code{type='location'}, either a single integer indicating the location in the data set to plot or a vector of length 2 providing the longitude and latitude of the new location.
#' @param loadings If \code{plot.param='lambda'}, an integer indicating which factor loading to plot.
#' @param factor If \code{type='factor'}, an integer or vector of integers specifying which factor(s) to plot.
#' @returns A plot of the selected parameter or variable.  
#' @author Candace Berrett and Adam Simpson
#' @examples
#' plot(out.sim)
#' @export 
plot.bstfa <- function(x, ..., plot.param='beta', type='map', location=1, loadings=NULL, factor=1){
  
  out <- x
  
  if(plot.param=="beta"){var.to.plot <- "slope"}
  if(plot.param=="mu"){var.to.plot <- "mean"}
  if(plot.param=="lambda"){var.to.plot <- "loading"}
  
  if(type=="map"){
    if(!plot.param %in% c("beta", "mu", "lambda")){stop("This type of plot requires a plot.param to be one of 'beta', 'mu', or 'lambda'.")}
     m <- map_spatial_param(out, parameter=var.to.plot, loadings=loadings, type='mean')
  }
  
  
  if(type=="basic"){
    if(plot.param != "beta" | plot.param != "mu" | plot.param != "lambda"){stop("This type of plot requires a plot.param to be one of 'beta', 'mu', or 'lambda'.")}
    m <- plot_spatial_param(out, parameter=var.to.plot, loadings=loadings, type='mean')
  }
  
  if(type=="location"){
    plot_location(out, location=location)
  }
  
  if(type=="annual"){
    plot_annual(out, location=location)
  }
  
  if(type=="factor"){
    plot_factor(out, factor=factor, together=TRUE)
  }
  
  invisible(NULL)
}


#' Print for 'bstfa' Objects
#'
#' This function prints the names and properties of a 'bstfa' object. 
#' @param x A bstfa object.
#' @param ... Additional arguments passed to or from other methods.
#' @returns A data frame containing the following elements:
#' \describe{
#'   \item{post.mean}{Posterior mean of parameters for each observed location.}
#'   \item{post.sd}{Posterior standard deviation of parameters for each observed location.}
#'   \item{post.pgt0}{Posterior probability that the parameter is greater than 0, if relevant.}
#' }
#' @author Candace Berrett and Adam Simpson
#' @examples
#' print(out.sim)
#' @export
print.bstfa <- function(x, ...){
  
  out <- x
  variable <- c("mu", "alpha.mu", "tau2.mu", "beta", "alpha.beta", 
                "tau2.beta", "xi", "alpha.xi", "tau2.xi", "F.tilde", 
                "alphaT", "Lambda.tilde", "alphaS", "tau2.lambda", 
                "sig2", "y.missing", "time.data", "setup.time", "model.matrices", 
                "factors.fixed", "iters", "y", "missing", "doy", "knots.spatial", 
                "knots.load", "draws")
  description <- c(paste("An mcmc object of size",  out$draws, "by",out$n.locs, "containing posterior draws for the mean of each location."),
  paste("An mcmc object of size", out$draws, "by",  out$n.spatial.bases, "containing posterior draws for the coefficients modeling the mean process."),
  paste("An mcmc object of size", out$draws, "by 1 containing the posterior draws for the variance of the mean process."),
  paste("An mcmc object of size", out$draws, "by", out$n.locs, "containing the posterior draws for the increase/decrease (slope) across time for each location."),
  paste("An mcmc object of size", out$draws, "by", out$n.spatial.bases, "containing posterior draws for the coefficients modeling the slope."),
  paste("An mcmc object of size", out$draws, "by 1 containing posterior draws of the variance of the slopes."),
  paste("An mcmc object of size", out$draws, "by", out$n.seasn.knots*out$n.locs, "containing posterior draws for the coefficients of the seasonal process."),
  paste("An mcmc object of size", out$draws, "by", out$n.spatial.bases*out$n.seasn.knots, "containing posterior draws for the coefficients modeling each coefficient of the seasonal process."),
  paste("An mcmc object of size", out$draws, "by", out$n.seasn.knots, "containing posterior draws of the variance of the coefficients of the seasonal process."),
  paste("An mcmc object of size", out$draws, "by", out$n.times*out$n.factors, "containing posterior draws of the residual factors."),
  paste("An mcmc object of size", out$draws, "by", out$n.factors*out$n.temp.bases, "containing posterior draws of the coefficients for the factor temporally-dependent process."),
  paste("An mcmc object of size", out$draws, "by", out$n.factors*out$n.locs, "containing posterior draws of the loadings for each location."),
  paste("An mcmc object of size", out$draws, "by", out$n.factors*out$n.load.bases, "containing posterior draws of the coefficients for the loadings spatial process."),
  paste("An mcmc object of size", out$draws, "by 1 containing posterior draws of the variance of the loadings spatial process."),
  paste("An mcmc object of size", out$draws, "by 1 containing posterior draws of the residual variance of the data."),
  ifelse(is.null(class(out$y.missing)), 'NULL', paste("A matrix of size", sum(out$missing), "by", out$draws, "containing posterior draws of the missing observations.")),
  paste("A data frame of size", out$iters, "by 6 containing the time it took to sample each parameter for every iteration."),
  paste("A difftime object containing the time the model setup took."),
  "A list containing the matrices used for each modeling process.",  
  paste("A vector of length", out$n.factors, "giving the location indices of the fixed loadings."),
  "A scalar returning the number of MCMC iterations.",
  paste("A vector of the", out$n.times*out$n.locs, "observations."),
  "A logical vector indicating whether that element's observation was missing or not.",
  paste("A numeric vector of length", out$n.times, "containing the day of year for each element in the original dates."),
  ifelse(out$spatial.style=="grid",  paste("A list of length", out$knot.levels, "containing the coordinates for all knots at each resolution."), "NA"),
  ifelse(out$load.style=="grid", paste("A list of length", out$knot.levels, "containing the coordinates for all knots at each resolution."), "NA"),
  "The number of saved MCMC iterations after removing the burn-in and thinning.")
  
  obj <- data.frame(variable=variable, description=description)
  
  for(jj in 1:nrow(obj)){
    cat(obj[jj,1], obj[jj,2], "\n")
  }
  
  return(invisible(x))
  
}



#' Coefficients for 'bstfa' Object
#'
#' This function provides the posterior mean values for main-effects coefficients. 
#' @param object A bstfa object.
#' @param ... Additional arguments passed to or from other methods.
#' @returns A numeric vector for BSTFA model coefficients.  Names of the values are labeled using the station name then the parameter.
#' @author Candace Berrett and Adam Simpson
#' @examples
#' coef(out.sim)
#' @export 
coef.bstfa <- function(object, ...){
  
  out <- object
  
  #Collect row names for data frame
  station.names <- row.names(out$coords)
  
  #posterior summaries of spatially-dependent mean
  if(out$mean){
    postm.mu <- apply(out$mu, 2, mean)
    df.names <- paste0(station.names, ".mu")
  }else{
    postm.mu <- NULL
    df.names <- NULL
  }
  
  #posterior summaries of spatially-dependent linear component
  if(out$linear){
    postm.beta <- apply(out$beta, 2, mean)
    df.names <- c(df.names, paste0(station.names, ".beta"))
  }else{
    postm.beta <- NULL
    df.names <- c(df.names, NULL)
  }
  
  #posterior summaries of seasonal spline coefficients
  if(out$seasonal){
    xi.nms <- expand.grid(paste0("xi.knot", 1:out$n.seasn.knots), station.names)
    xi.nms <- paste0(xi.nms$Var2, ".",  xi.nms$Var1)
    postm.xi <- apply(out$xi, 2, mean)
    df.names <- c(df.names, xi.nms)
  }else{
    postm.xi <- NULL
  }
  
  coefficients <- c(postm.mu, postm.beta, postm.xi)
  names(coefficients) <- df.names
  
  return(coefficients)
  
}
