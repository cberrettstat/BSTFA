##### BSTFA FUNCTION - FA Reduction built in #####

#' Reduced BSTFA function
#'
#' This function uses MCMC to draw from posterior distributions of a Bayesian spatio-temporal factor analysis model.  All spatial processes use one of Fourier, thin plate spline, or multiresolution basis functions.  The temporally-dependent factors use Fourier bases.  The default values are chosen to work well for many data sets.  Thus, it is possible to use this function using only three arguments: \code{ymat}, \code{dates}, and \code{coords}.  The default number of MCMC iterations is 10000 (saving 5000); however, depending on the number of observations and processes modeled, it may need more draws than this to ensure the posterior draws are representative of the entire posterior distribution space.
#' @param ymat Data matrix of size \code{n.times} by \code{n.locs}. Any missing data should be marked by \code{NA}.  The model works best if the data are zero-centered for each location.
#' @param dates \code{n.times} length vector of class \code{'Date'} corresponding to each date of the observed data.  For now, the dates should be regularly spaced (e.g., daily).
#' @param coords \code{n.locs} by \code{2} matrix or data frame of coordinates for the locations of the observed data. If using longitude and latitude, longitude is assumed to be the first coordinate.
#' @param iters Number of MCMC iterations to draw.  Default value is \code{10000}.  Function only saves \code{(iters-burn)/thin} drawn values.
#' @param n.times Number of observations for each location. Default is \code{nrow(ymat)}.
#' @param n.locs Number of observed locations.  Default is \code{ncol(ymat)}.
#' @param x Optional \code{n.locs} by \code{p} matrix of covariates for each location.  If there are no covariates, set to \code{NULL} (default).
#' @param mean Logical scalar.  If \code{TRUE}, the model will fit a spatially-dependent mean for each location.  Otherwise, the model will assume the means are zero at each location (default).
#' @param linear Logical scalar.  If \code{TRUE} (default), the model will fit a spatially-dependent linear increase/decrease (or slope) in time. Otherwise, the model will assume a zero change in slope across time.
#' @param seasonal Logical scalar. If \code{TRUE} (default), the model will use circular b-splines to model a spatially-dependent annual process.  Otherwise, the model will assume there is no seasonal (annual) process.
#' @param factors Logical scalar. If \code{TRUE} (default), the model will fit a spatio-temporal factor analysis model with temporally-dependent factors and spatially-dependent loadings.
#' @param n.seasn.knots Numeric scalar indicating the number of knots to use for the seasonal basis components. The default value is \code{min(7, ceiling(length(unique(yday(dates)))/3))}, where 7 will capture approximately 2 peaks during the year.
#' @param spatial.style Character scalar indicating the style of bases to use for the linear and seasonal components.  Style options are \code{'fourier'} (default), \code{'tps'} for thin plate splines, and \code{'grid'} for multiresolution bisquare bases using knots from a grid across the space.
#' @param n.spatial.bases Numeric scalar indicating the number of spatial bases to use when \code{spatial.style} is either \code{'fourier'} or \code{'tps'}. Default value is \code{min(16, ceiling(n.locs/3))}. When \code{spatial.style} is \code{'fourier'}, this value must be an even square number.
#' @param knot.levels Numeric scalar indicating the number of resolutions to use for when \code{spatial.style='grid'} and/or \code{load.style='grid'}.  Default is 2.
#' @param max.knot.dist Numeric scalar indicating the maximum distance at which a basis value is greater than zero when \code{spatial.style='grid'} and/or \code{load.style='grid'}.  Default value is \code{mean(dist(coords))}.
#' @param premade.knots Optional list of length \code{knot.levels} with each list element containing a matrix of longitude-latitude coordinates of the knots to use for each resolution when \code{spatial.style='grid'} and/or \code{load.style='grid'}.  Otherwise, when \code{premade.knots = NULL} (default), the knots are determined by using the standard multiresolution grids across the space.
#' @param plot.knots Logical scalar indicating whether to plot the knots used when \code{spatial.style='grid'} and/or \code{load.style='grid'}. Default is \code{FALSE}.
#' @param n.factors Numeric scalar indicating how many factors to use in the model.  Default is \code{min(4,ceiling(n.locs/20))}.
#' @param factors.fixed Numeric vector of length \code{n.factors} indicating the locations to use for the fixed loadings.  This is needed for model identifiability.  If \code{factors.fixed=NULL} (default), the code will select locations with less than 20% missing data and that are far apart in the space.
#' @param plot.factors Logical scalar indicating whether to plot the fixed factor locations.  Default is \code{FALSE}.
#' @param load.style Character scalar indicating the style of spatial bases to use for the spatially-dependent loadings. Options are \code{'fourier'} (default) for the Fourier bases, \code{'tps'} for thin plate splines, and \code{'grid'} for multiresolution bases.  This can be the same as or different than \code{spatial.style}.
#' @param n.load.bases Numeric scalar indicating the number of bases to use for the spatially-dependent loadings when \code{load.style} is either \code{'fouier'} or \code{'tps'}.  This can be the same as or different than  \code{n.spatial.bases}.  Default is \code{4}. When \code{load.style='fourier'}, this value must be an even square number.
#' @param freq.lon Numeric scalar used for \code{spatial.style} or \code{load.style} equal to \code{'fourier'} or \code{'eigen'}.  For \code{'fourier'}, this is the frequency used for the first column of \code{coords} (assumed to be longitude) for the Fourier bases. For \code{'eigen'}, this is the range parameter of the exponential spatial correlation matrix used to create the eigenvectors. Default value is \code{2*diff(range(coords[,1]))}.
#' @param freq.lat Numeric scalar used for \code{spatial.style} or \code{load.style} equal to \code{'fourier'}.  This is the frequency to use for the second column of \code{coords} (assumed to be latitude) for the Fourier bases. Default value is \code{2*diff(range(coords[,2]))}.
#' @param n.temp.bases Numeric scalar indicating the number of Fourier bases to use for the temporally-dependent factors. The default value is 10% of \code{n.times}.
#' @param freq.temp Numeric scalar indicating the frequency to use for the Fourier bases of the temporally-dependent factors.  The default value is \code{n.times}.
#' @param alpha.prec Numeric scalar indicating the prior precision for all model process coefficients. Default value is \code{1/100000}.
#' @param tau2.gamma Numeric scalar indicating the prior shape for the precision of the model coefficients.  Default value is \code{2}.
#' @param tau2.phi Numeric scalar indicating the prior rate for the precision of the model coefficients.  Default value is \code{1e-07}.
#' @param sig2.gamma Numeric scalar indicating the prior shape for the residual precision.  Default value is \code{2}.
#' @param sig2.phi Numeric scalar indicating the prior rate for the residual precision. Default value is \code{1e-05}.
#' @param sig2 Numeric scalar indicating the starting value for the residual variance. If \code{NULL} (default), the function will select a reasonable starting value.
#' @param beta Numeric vector of length \code{n.locs + p} indicating starting values for the slopes.  If \code{NULL} (default), the function will select reasonable starting values.
#' @param xi Numeric vector of length \code{(n.locs + p)*n.seasn.knots} indicating starting values for the coefficients of the seasonal component. If \code{NULL} (default), the function will select reasonable starting values.
#' @param Fmat Numeric matrix of size \code{n.times} by \code{n.factors} indicating starting values for the factors.  Default value is to start all factor values at 0.
#' @param Lambda Numeric matrix of size \code{n.locs} by \code{n.factors} indicating starting values for the loadings.  Default value is to start all loadings at 0.
#' @param thin Numeric scalar indicating how many MCMC iterations to thin by.  Default value is 1, indicating no thinning.
#' @param burn Numeric scalar indicating how many MCMC iterations to burn before saving.  Default value is one-half of \code{iters}.
#' @param verbose Logical scalar indicating whether or not to print the status of the MCMC process.  If \code{TRUE} (default), the function will print every time an additional 10% of the MCMC process is completed.
#' @param save.missing Logical scalar indicating whether or not to save the MCMC draws for the missing observations.  If \code{TRUE} (default), the function will save an additional MCMC object containing the MCMC draws for each missing observation.  Use \code{FALSE} to save file space and memory.
#' @param save.time Logical scalar indicating whether to save the computation time for each MCMC iteration.  Default value is \code{FALSE}.  When \code{FALSE}, the function \code{compute_summary()} will not be useful.
#' @param marginalize Logical scalar indicating whether to sample hyper-coefficients by marginalizing out the corresponding parameters.  Default value is \code{TRUE}.
#' @importFrom matrixcalc vec
#' @importFrom mgcv cSplineDes
#' @importFrom coda as.mcmc
#' @importFrom MASS mvrnorm
#' @importFrom npreg basis.tps
#' @importFrom lubridate yday
#' @importFrom utils combn
#' @import matrixcalc
#' @import npreg
#' @import stats
#' @import graphics
#' @returns A list containing the following elements (any elements that are the same as in the function input are removed here for brevity):
#' \describe{
#'   \item{mu}{An mcmc object of size \code{draws} by \code{n.locs} containing posterior draws for the mean of each location.  If \code{mean=FALSE} (default), the values will all be zero.}
#'   \item{alpha.mu}{An mcmc object of size \code{draws} by \code{n.spatial.bases + p} containing posterior draws for the coefficients modeling the mean process.  If \code{mean=FALSE} (default), the values will all be zero.}
#'   \item{tau2.mu}{An mcmc object of size \code{draws} by \code{1} containing the posterior draws for the variance of the mean process.  If \code{mean=FALSE} (default), the values will all be zero.}
#'   \item{beta}{An mcmc object of size \code{draws} by \code{n.locs} containing the posterior draws for the increase/decrease (slope) across time for each location.}
#'   \item{alpha.beta}{An mcmc object of size \code{draws} by \code{n.spatial.bases + p} containing posterior draws for the coefficients modeling the slope.}
#'   \item{tau2.beta}{An mcmc object of size \code{draws} by \code{1} containing posterior draws of the variance of the slopes.}
#'   \item{xi}{An mcmc object of size \code{draws} by \code{n.seasn.knots*n.locs} containing posterior draws for the coefficients of the seasonal process.}
#'   \item{alpha.xi}{An mcmc object of size \code{draws} by \code{(n.spatial.bases + p)*n.seasn.knots} containing posterior draws for the coefficients modeling each coefficient of the seasonal process.}
#'   \item{tau2.xi}{An mcmc object of size \code{draws} by \code{n.seasn.knots} containing posterior draws of the variance of the coefficients of the seasonal process.}
#'   \item{F.tilde}{An mcmc object of size \code{draws} by \code{n.times*n.factors} containing posterior draws of the residual factors.}
#'   \item{alphaT}{An mcmc object of size \code{draws} by \code{n.factors*n.temp.bases} containing posterior draws of the coefficients for the factor temporally-dependent process.}
#'   \item{Lambda.tilde}{An mcmc object of size \code{draws} by \code{n.factors*n.locs} containing posterior draws of the loadings for each location.}
#'   \item{alphaS}{An mcmc object of size \code{draws} by \code{n.factors*n.load.bases} containing posterior draws of the coefficients for the loadings spatial process.}
#'   \item{tau2.lambda}{An mcmc object of size \code{draws} by \code{1} indicating the residual variance of the loadings spatial process.}
#'   \item{sig2}{An mcmc object of size \code{draws} by \code{1} containing posterior draws of the residual variance of the data.}
#'   \item{y.missing}{If \code{save.missing=TRUE}, a matrix of size \code{sum(missing)} by \code{draws} containing posterior draws of the missing observations.  Otherwise, the object is \code{NULL}. }
#'   \item{time.data}{A data frame of size \code{iters} by \code{6} containing the time it took to sample each parameter for every iteration.}
#'   \item{setup.time}{An object containing the time the model setup took.}
#'   \item{model.matrices}{A list containing the matrices used for each modeling process. \code{newS} is the matrix of spatial basis coefficients for the mean, linear, and seasonal process coefficients.  \code{linear.Tsub} is the matrix used to enforce a linear increase/increase (slope) across time. \code{seasonal.bs.basis} is the matrix containing the circular b-splines of the seasonal process.  \code{confoundingPmat.prime} is the matrix that enforces orthogonality of the factors from the mean, linear, and seasonal processes.  \code{QT} contains the Fourier bases used to model the temporal factors.  \code{QS} contains the bases used to model the spatial loadings.}
#'   \item{factors.fixed}{A vector of length \code{n.factors} giving the location indices of the fixed loadings.}
#'   \item{iters}{A scalar returning the number of MCMC iterations.}
#'   \item{y}{An \code{n.times*n.locs} vector of the observations.}
#'   \item{missing}{A logical vector indicating whether that element's observation was missing or not.}
#'   \item{doy}{A numeric vector of length \code{n.times} containing the day of year for each element in the original \code{dates}.}
#'   \item{knots.spatial}{For \code{spatial.style='grid'}, a list of length \code{knot.levels} containing the coordinates for all knots at each resolution.}
#'   \item{knots.load}{For \code{load.style='grid'}, a list of length \code{knot.levels} containing the coordinates for all knots at each resolution.}
#'   \item{draws}{The number of saved MCMC iterations after removing the burn-in and thinning.}
#' }
#' @author Adam Simpson and Candace Berrett
#' @examples
#' #Very small example to illustrate use and ensure functionality
#' data(utahDataList)
#' attach(utahDataList)
#'
#' #identify locations with very little missing data just for this example
#' low.miss <- which(apply(is.na(TemperatureVals), 2, mean)<.02)
#'
#' out <- BSTFA(ymat=TemperatureVals[1:50,low.miss],
#'   dates=Dates[1:50],
#'   coords=Coords[low.miss,],
#'   n.factors=2,
#'   iters=10)
#'
#' # More full example:
#' # out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords)
#' @export BSTFA
BSTFA <- function(ymat, dates, coords,
                 iters=10000, n.times=nrow(ymat), n.locs=ncol(ymat), x=NULL,
                 mean=FALSE, linear=TRUE, seasonal=TRUE, factors=TRUE,
                 n.seasn.knots=min(7, ceiling(length(unique(yday(dates)))/3)),
                 spatial.style='eigen',
                 n.spatial.bases=min(10, ceiling(n.locs/3)),
                 knot.levels=2, max.knot.dist=mean(dist(coords)), premade.knots=NULL, plot.knots=FALSE,
                 n.factors=min(4,ceiling(n.locs/20)), factors.fixed=NULL, plot.factors=FALSE,
                 load.style='eigen',
                 n.load.bases=4,
                 freq.lon=2*diff(range(coords[,1])),
                 freq.lat=2*diff(range(coords[,2])),
                 n.temp.bases=ifelse(floor(n.times*0.10)%%2==1, floor(n.times*0.10)-1, floor(n.times*0.10)),
                 freq.temp=n.times,
                 alpha.prec=1/100000, tau2.gamma=2, tau2.phi=0.0000001, sig2.gamma=2, sig2.phi=1e-5,
                 sig2=NULL, beta=NULL, xi=NULL,
                 Fmat=matrix(0,nrow=n.times,ncol=n.factors), Lambda=matrix(0,nrow=n.locs, n.factors),
                 thin=1, burn=iters*0.5, verbose=TRUE, save.missing=TRUE, save.time=FALSE, 
                 marginalize=TRUE) {

  
  start <- Sys.time()

  #Basic checks on inputs
  if(n.spatial.bases > n.locs){stop("n.spatial.bases must be less than n.locs")}
  if(n.load.bases > n.locs){stop("n.load.bases must be less than n.locs")}
  if(n.temp.bases > n.times){stop("n.temp.bases must be less than n.times")}
  
  if(!is.matrix(ymat)){ymat <- as.matrix(ymat)}
  if(!is.null(factors.fixed)){n.factors<-length(factors.fixed)}

  ### Prepare missing data
  # Make missing values 0 for now, but they will be estimated differently
  y <- c(ymat)
  missing = ifelse(is.na(y), TRUE, FALSE)
  prop.missing = apply(ymat, 2, function(x) sum(is.na(x)) / n.times)
  missind <- which(missing)
  notmissind <- which(!missing)
  y[missing] = 0
  if (is.null(sig2)) sig2 = var(y)/10

  if(save.missing==T & sum(missing)!=0){
    y.save <- matrix(0, nrow=sum(missing), ncol=floor((iters-burn)/thin))
  }else{
    y.save <- NULL
  }

  ### Create doy
  doy <- yday(dates)

  ### Change x to matrix if not null
  if (!is.null(x)) x <- as.matrix(x)

  ### Change coordinates to a matrix if not
  if(!is.matrix(coords)){coords <- as.matrix(coords)}
  ### Create newS
  knots.vec.save.spatial=NULL
  knots.vec.save.load=NULL
  if (spatial.style=='grid') {
    if (is.null(premade.knots)) {
      n.spatial.bases=0
      for (i in 1:(knot.levels)) {
        n.spatial.bases <- n.spatial.bases + 4^(i)
      }
    }
    else {
      n.spatial.bases = nrow(matrix(unlist(premade.knots),ncol=2))
    }
    ### using function makeNewS - uses bisquare distance
    newS.output = makeNewS(coords=coords,n.locations=n.locs,knot.levels=knot.levels,
                           max.knot.dist=max.knot.dist, x=x,
                           plot.knots=plot.knots,
                           premade.knots=premade.knots)

    newS = newS.output[[1]]
    knots.vec.save.spatial = newS.output[[2]]
  }
  if (spatial.style=='fourier') {
    if (sqrt(n.spatial.bases)%%1 != 0 | n.spatial.bases%%2 != 0) {
      n.spatial.bases=floor(sqrt(n.spatial.bases))^2
      if(n.spatial.bases%%2 != 0){n.spatial.bases <- (floor(sqrt(n.spatial.bases))+1)^2}
      message(paste("n.spatial.bases must be an even square number; changed value to", n.spatial.bases))
    }
    ### Original Fourier Method
    m.fft.lon <- sapply(1:(sqrt(n.spatial.bases)/2), function(k) {
      sin_term <- sin(2 * pi * k * (coords[,1])/freq.lon)
      cos_term <- cos(2 * pi * k * (coords[,1])/freq.lon)
      cbind(sin_term, cos_term)
    })
    m.fft.lat <- sapply(1:(sqrt(n.spatial.bases)/2), function(k) {
      sin_term <- sin(2 * pi * k * (coords[,2])/freq.lat)
      cos_term <- cos(2 * pi * k * (coords[,2])/freq.lat)
      cbind(sin_term, cos_term)
    })
    Slon <- cbind(m.fft.lon[1:n.locs,], m.fft.lon[(n.locs+1):(2*n.locs),])
    Slat <- cbind(m.fft.lat[1:n.locs,], m.fft.lat[(n.locs+1):(2*n.locs),])
    newS <- matrix(NA, nrow=n.locs, ncol=n.spatial.bases)
    col_idx <- 1
    for (thisi in 1:ncol(Slon)) {
      for (thisj in 1:ncol(Slat)) {
        newS[, col_idx] <- Slon[, thisi] * Slat[, thisj]
        col_idx <- col_idx + 1
      }
    }
    
    if(!is.null(x)){
      newS <- cbind(newS, x)
    }
    if (qr(newS)$rank != ncol(newS)) {
      stop("Collinearity in bases for spatial coefficients; adjust Fourier frequencies.")
    }
    
  }
  
  if(spatial.style=='eigen'){
    distmat <- as.matrix(dist(coords))
    cormat <- exp(-distmat/freq.lon)
    eigs <- eigen(cormat)
    newS <- eigs$vectors[,1:n.spatial.bases]
    if (!is.null(x)){
      newS <- cbind(newS, x)
      A.prec <- diag(c(.001*1/eigs$values[1:n.spatial.bases], rep(alpha.prec, dim(x)[2])))
    }else{
      A.prec <- .001*solve(diag(eigs$values[1:n.spatial.bases]))
    }
  }
  
  if (spatial.style=='tps') {
    dd = floor(sqrt(n.spatial.bases))
    xxx = yyy = seq(1/(dd*2), (dd*2-1)/(dd*2), by=1/dd)
    knots.spatial = expand.grid(xxx,yyy)
    range_long = max(coords[,1]) - min(coords[,1])
    range_lat = max(coords[,2]) - min(coords[,2])
    knots.spatial[,1] = (knots.spatial[,1] * range_long) + min(coords[,1])
    knots.spatial[,2] = (knots.spatial[,2] * range_lat) + min(coords[,2])
    knots.vec.save.spatial = knots.spatial
    newS = npreg::basis.tps(coords,
                     knots=knots.spatial,
                     rk=FALSE)[,-(1:2)]
    n.spatial.bases = dd^2
    if(!is.null(x)){newS <- cbind(newS,x)}
  }
  
  if(spatial.style!='eigen'){
    A.prec = diag(alpha.prec, dim(newS)[2])
  }

  model.matrices <- list()
  model.matrices$newS <- newS
  model.matrices$A.prec <- A.prec
  
  ### Set up mean component
  if(mean==TRUE){
    Jfull = kronecker(Matrix::Diagonal(n=n.locs), rep(1, n.times))
    Jfullmiss = Jfull[notmissind,,drop=FALSE]
    tJ <- Matrix::t(Jfullmiss)
    tJJ <- crossprod(Jfullmiss)
    #ItJJ <- methods::as(kronecker(diag(1,n.locs), t(rep(1,n.times))%*%rep(1,n.times)), "sparseMatrix")
    #ItJ <- methods::as(kronecker(diag(1,n.locs), t(rep(1,n.times))), "sparseMatrix")
    mu.var <- solve(tJJ) #solve(ItJJ)
    mu.mean <- mu.var%*%tJ%*%y[-missind] #mu.var%*%ItJ%*%y
    mu <-  my_mvrnorm(mu.mean, 0.001*mu.var)
    mu <- as.matrix(mu)
    Jfullmu.long <- Jfull%*%mu
    rm(list=c("mu.mean", "mu.var"))
    alpha.mu=rep(0, dim(newS)[2])
    tau2.mu = as.numeric(var(mu))
  } else {
    mu <- rep(0, n.locs)
    Jfullmu.long <- rep(0, n.times*n.locs)
  }
  mu.save <- matrix(0, nrow=n.locs, ncol=floor((iters-burn)/thin))
  alpha.mu.save <- matrix(0, nrow=dim(newS)[2], ncol=floor((iters-burn)/thin))
  tau2.mu.save <- matrix(0,nrow=1,ncol=floor((iters-burn)/thin))


  ### Set up linear component
  if (linear == TRUE) {
    Tsub <- -(n.times/2-0.5):(n.times/2-0.5)
    Tfull <- kronecker(Matrix::Diagonal(n=n.locs), Tsub)
    Tfullmiss <- Tfull[notmissind,,drop=FALSE]
    tT <- Matrix::t(Tfullmiss)
    tTT <- crossprod(Tfullmiss) #tT%*%Tfullmiss
    #ItTT <- methods::as(kronecker(Matrix::Diagonal(n=n.locs), t(Tsub)%*%Tsub), "sparseMatrix")
    #ItT <- methods::as(kronecker(Matrix::Diagonal(n=n.locs), t(Tsub)), "sparseMatrix")
    if(is.null(beta)==T){
      beta.var <- solve(tTT) #solve(ItTT)
      beta.mean <- beta.var%*%tT%*%(y[notmissind] - Jfullmu.long[notmissind]) #ItT%*%y #starting values for beta
      beta <- my_mvrnorm(beta.mean, 0.001*beta.var)
      beta <- as.matrix(beta)
      #beta <- beta + rnorm(length(beta.mean), 0, sd(beta.mean))
      rm(list=c("beta.mean", "beta.var"))
    }
    Tfullbeta.long <- Tfull%*%beta
    model.matrices$linear.Tsub <- Tsub
    alpha.beta <- rep(0, dim(newS)[2])
    tau2.beta <- as.numeric(var(beta))
 } else {
    beta <- rep(0, n.locs)
    Tfullbeta.long <- rep(0, n.times*n.locs)
  }
  beta.save <- matrix(0, nrow=n.locs, ncol=floor((iters-burn)/thin))
  alpha.beta.save <- matrix(0, nrow=dim(newS)[2], ncol=floor((iters-burn)/thin))
  tau2.beta.save <- matrix(0,nrow=1,ncol=floor((iters-burn)/thin))

  ######### If eigen style, can use variog and variofit to estimate phi get a "better" newS?
  

  ### Set up seasonal component
  model.matrices$seasonal.bs.basis <- matrix(0,nrow=n.times,ncol=n.seasn.knots)
  if(seasonal == TRUE) {
    newS.xi <- methods::as(kronecker(newS, diag(n.seasn.knots)), "sparseMatrix")
    # newS.xi <- kronecker(newS,diag(n.seasn.knots))
    knots <- seq(1, 366, length=n.seasn.knots+1)
    bs.basis <- mgcv::cSplineDes(doy, knots)
    Bfull <- kronecker(Matrix::Diagonal(n=n.locs), bs.basis)
    Bfullmiss <- Bfull[notmissind,,drop=FALSE]
    tB <- Matrix::t(Bfullmiss)
    tBB <- crossprod(Bfullmiss) #tB%*%Bfullmiss
    #ItBB <- methods::as(kronecker(Matrix::Diagonal(n=n.locs), t(bs.basis)%*%bs.basis), "sparseMatrix")
    #ItB <- methods::as(kronecker(Matrix::Diagonal(n=n.locs), t(bs.basis)), "sparseMatrix")
    if (is.null(xi)) {
      xi.var <- solve(tBB) #solve(ItBB)
      xi.mean <- xi.var%*%tB%*%(y[notmissind]-Jfullmu.long[notmissind]-Tfullbeta.long[notmissind]) #ItB%*%(y - Tfullbeta.long)
      xi <- my_mvrnorm(xi.mean, 0.001*xi.var) #+ rnorm(length(xi.mean), 0, sd(xi.mean)) #starting values for xi
      rm(list=c("xi.var", "xi.mean"))
    }
    Bfullxi.long <- Bfull%*%xi
    model.matrices$seasonal.bs.basis <- bs.basis
    alpha.xi <- rep(0, dim(newS.xi)[2])
    tau2.xi <- apply(matrix(xi, nrow=n.seasn.knots, ncol=n.locs, byrow=F), 1, var)/sqrt(n.locs) #rep(0.01, n.seasn.knots) #as.numeric(var(xi))
  } else {
    xi <- rep(0, n.locs*n.seasn.knots)
    Bfullxi.long <- rep(0, n.locs*n.times)
  }
  xi.save <- matrix(0, nrow=n.locs*n.seasn.knots, ncol=floor((iters-burn)/thin))
  alpha.xi.save <- matrix(0, nrow=n.seasn.knots*dim(newS)[2], ncol=floor((iters-burn)/thin))
  tau2.xi.save <- matrix(0, nrow=n.seasn.knots, ncol=floor((iters-burn)/thin))


  ### Deal with confounding
  if (mean | linear | seasonal) {
    Cmat <- NULL
    if (mean) Cmat <- cbind(Cmat, rep(1,n.times))
    if (linear) Cmat <- cbind(Cmat, Tsub)
    if (seasonal) Cmat <- cbind(Cmat, bs.basis)
    tCC <- crossprod(Cmat) 
    tCC <- (t(tCC) + tCC)/2
    if (mean) {
      Pmat <- Cmat%*%MASS::ginv(tCC)%*%t(Cmat)
    } else {
      Pmat <- Cmat%*%solve(tCC)%*%t(Cmat)
    }
    Pmat.prime <- diag(1, n.times) - Pmat
  } else {
    Pmat.prime = diag(1, n.times)
  }
  model.matrices$confoundingPmat.prime = Pmat.prime


  ### Set up Factor Analysis
  
  if (factors) {
    ### Set up temporal FA

    ### Eigen Method
    # if (is.null(phi.T)) phi.T = n.times/2
    # distT <- as.matrix(dist(1:n.times))
    # corT <- exp(-distT/phi.T)
    # QT <- eigen(corT)$vectors[,1:n.temp.bases]
    # bigQT <- kronecker(QT, diag(1,n.factors))

    ### Fourier Method
    # Create Fourier basis functions
    if (n.temp.bases%%2 == 1) n.temp.bases=n.temp.bases+1
    m.fft <- sapply(1:(n.temp.bases/2), function(k) {
      sin_term <- sin(2 * pi * k * (1:n.times)/freq.temp)
      cos_term <- cos(2 * pi * k * (1:n.times)/freq.temp)
      cbind(sin_term, cos_term)
    })
    QT <- cbind(m.fft[1:n.times,], m.fft[(n.times+1):(2*n.times),])
    model.matrices$QT = QT

    ##################################

    ### Set up spatial FA
    tau2.lambda=0.01

    ### Bisquare Method
    if (load.style == 'grid') {
      if (is.null(premade.knots)) {
        n.load.bases=0
        for (i in 1:(knot.levels)) {
          n.load.bases <- n.load.bases + 4^(i)
        }
      }
      else {
        n.load.bases = nrow(matrix(unlist(premade.knots),ncol=2))
      }
      if (spatial.style=='grid') {
        QS = newS[,1:n.load.bases]
        knots.vec.save.load = newS.output[[2]]
      }
      else {
        newS.output = makeNewS(coords=coords,n.locations=n.locs,knot.levels=knot.levels,
                               max.knot.dist=max.knot.dist, x=NULL,
                               plot.knots=plot.knots,
                               premade.knots=premade.knots)
        QS = newS.output[[1]]
        knots.vec.save.load = newS.output[[2]]
      }
    }

    ### Fourier method
    if (load.style == 'fourier') {
      if (sqrt(n.load.bases)%%1 != 0 | n.load.bases%%2 != 0) {
        n.load.bases=floor(sqrt(n.load.bases))^2
        if(n.load.bases%%2 != 0){n.load.bases <- (floor(sqrt(n.load.bases))+1)^2}
        message(paste("n.load.bases must be an even square number; changed value to", n.load.bases))
      }
      ### Original Fourier Method
      m.fft.lon <- sapply(1:(sqrt(n.load.bases)/2), function(k) {
        sin_term <- sin(2 * pi * k * (coords[,1])/freq.lon)
        cos_term <- cos(2 * pi * k * (coords[,1])/freq.lon)
        cbind(sin_term, cos_term)
      })
      m.fft.lat <- sapply(1:(sqrt(n.load.bases)/2), function(k) {
        sin_term <- sin(2 * pi * k * (coords[,2])/freq.lat)
        cos_term <- cos(2 * pi * k * (coords[,2])/freq.lat)
        cbind(sin_term, cos_term)
      })
      QSlon <- cbind(m.fft.lon[1:n.locs,], m.fft.lon[(n.locs+1):(2*n.locs),])
      QSlat <- cbind(m.fft.lat[1:n.locs,], m.fft.lat[(n.locs+1):(2*n.locs),])
      QS <- matrix(NA, nrow=n.locs, ncol=n.load.bases)
      col_idx <- 1
      for (thisi in 1:ncol(QSlon)) {
        for (thisj in 1:ncol(QSlat)) {
          QS[, col_idx] <- QSlon[, thisi] * QSlat[, thisj]
          col_idx <- col_idx + 1
        }
      }
      
      if (qr(QS)$rank != ncol(QS)) {
        stop("Collinearity in bases for spatial loadings; adjust Fourier frequencies.")
      }
      
    }
    
    if(load.style=='eigen'){
      distmat <- as.matrix(dist(coords))
      cormat <- exp(-distmat/freq.lon)
      eigs <- eigen(cormat)
      QS <- eigs$vectors[,1:n.load.bases]
      A.lambda.prec <- solve(diag(eigs$values[1:n.load.bases]))
    }
    
    if (load.style=='tps') {
      dd = floor(sqrt(n.load.bases))
      xxx = yyy = seq(1/(dd*2), (dd*2-1)/(dd*2), by=1/dd)
      knots.load = expand.grid(xxx,yyy)
      range_long = max(coords[,1]) - min(coords[,1])
      range_lat = max(coords[,2]) - min(coords[,2])
      knots.load[,1] = (knots.load[,1] * range_long) + min(coords[,1])
      knots.load[,2] = (knots.load[,2] * range_lat) + min(coords[,2])
      knots.vec.save.load = knots.load
      QS = npreg::basis.tps(coords,
                     knots=knots.load,
                     rk=FALSE)[,-(1:2)]
      n.load.bases = dd^2
    }
    
    if(load.style!='eigen'){
      A.lambda.prec = diag(alpha.prec, dim(QS)[2])
    }
    
    model.matrices$A.lambda.prec <- A.lambda.prec

    ### Gaussian Method

    # gaussian_basis <- function(x, y, mu_x, mu_y, sigma) {
    #   return(exp(-((x - mu_x)^2 + (y - mu_y)^2) / (2 * sigma^2)))
    # }
    # # Parameters for multiple Gaussians
    # ### FIX the "by=1", and the "sigma=1" as parameters in the function
    # bumps = expand.grid(seq(min(coords[,1]), max(coords[,1]), length.out=5),
    #                     seq(min(coords[,2]), max(coords[,2]), length.out=5))
    # sigma=1
    # n.load.bases = 5*5
    # QS = matrix(0,nrow=nrow(coords),ncol=nrow(bumps))
    # # Sum the Gaussian basis functions
    # for (i in 1:nrow(bumps)) {
    #   QS[,i] = gaussian_basis(coords[,1], coords[,2], bumps[i,1], bumps[i,2], sigma)
    # }
    # model.matrices$QS = QS

    model.matrices$QS = QS


    ##################################

    ### Establish fixed factor locations
    if (is.null(factors.fixed)) {
      if(n.factors==1){
        factors.fixed <- sample(which(prop.missing==min(prop.missing, na.rm=T)), size=1)
      }else{
      distmat <- as.matrix(dist(coords))
      far = FALSE
      d = c()
      while (!far) {
        p = sample(which(prop.missing<0.2), size=n.factors, replace=FALSE)
        combos = combn(p,2)
        for (i in 1:ncol(combos)) {
          d[i] = distmat[combos[1,i], combos[2,i]]
        }
        far = ifelse(min(d) < (max(distmat) / n.factors), FALSE, TRUE)
      }
      factors.fixed = p
      }
    }
    Lambda.tilde = Lambda
    Lambda.tilde[factors.fixed,] = diag(n.factors)

    # if (plot.factors) {
    #   plot(coords, xlab='Longitude', ylab='Latitude', main='Fixed Factor Locations')
    #   points(coords[factors.fixed,], col='red', cex=2, pch=19)
    # }
  }

  delayFA = min(floor(burn/2), 500)

  alphaT.save <- matrix(0, nrow=n.factors*n.temp.bases, ncol=floor((iters-burn)/thin))
  F.tilde.save <- matrix(0, nrow=n.factors*n.times, ncol=floor((iters-burn)/thin))
  Lambda.tilde.save <- matrix(0, nrow=n.factors*n.locs, ncol=floor((iters-burn)/thin))
  alphaS.save <- matrix(0, nrow=n.factors*n.load.bases, ncol=floor((iters-burn)/thin))
  tau2.lambda.save <- matrix(0, nrow=1, ncol=floor((iters-burn)/thin))
  alphaT <- rep(0, n.factors*n.temp.bases)
  alphaS <- rep(0, n.factors*n.load.bases)
  FLambda.long <- c(Fmat%*%t(Lambda))
  
  ### Set up variance component
  sig2.save <- matrix(0, nrow=1, ncol=floor((iters-burn)/thin))

  ### Useful one-time calculations
  
  if (seasonal){
    StSI <- methods::as(kronecker(t(newS)%*%newS, Matrix::Diagonal(n=n.seasn.knots)), "sparseMatrix")
  }
  if (factors) {
    PQT <- Pmat.prime%*%QT
    PQTtPQT = t(PQT)%*%PQT
    QsI <- methods::as(kronecker(QS, diag(1, n.factors)), "sparseMatrix")
    #QsI <- kronecker(QS, diag(1,n.factors))
    QstQsI <- methods::as(kronecker(t(QS)%*%QS, diag(1, n.factors)), "sparseMatrix")
    #QstQsI <- kronecker(t(QS)%*%QS, diag(1, n.factors))
  }
  if(linear){
    TtT <- crossprod(tT)
  }
  if(seasonal){
    BtB <- crossprod(tB)
  }
  
  ### Set up effective sample size calculations
  eSS.check=1000
  eSS.converged=100

  ### Set up time.data
  time.data = matrix(0, nrow=floor(iters/thin), ncol=5)
  time.data = as.data.frame(time.data)
  colnames(time.data) = c('beta', 'xi', 'F.tilde', 'Lambda.tilde', 'sigma2')
  end <- Sys.time()
  setup.time = end-start

  if (verbose) cat(paste("Setup complete! Time taken: ", round(setup.time/60,2), " minutes. \n", sep=""))
  if (verbose) cat(paste("Starting MCMC, ", iters, " iterations. \n", sep=""))

  ### MCMC ###
  start.time = proc.time()
  for(i in 1:iters){

    ### Sample values of mu
    if (mean) {
      temp <- y - Tfullbeta.long - Bfullxi.long - FLambda.long
      mu.var <- solve((1/sig2)*tJJ + (1/tau2.mu)*Matrix::Diagonal(n=n.locs)) #ItJJ + (1/tau2.mu)*Matrix::Diagonal(n=n.locs))
      mu.mean <- mu.var%*%((1/sig2)*tJ%*%temp[notmissind] + (1/tau2.mu)*newS%*%alpha.mu) #ItJ%*%temp + (1/tau2.mu)*newS%*%alpha.mu)
      mu <- as.vector(MASS::mvrnorm(1,mu.mean,mu.var))
      Jfullmu.long <- Jfull%*%mu
      rm(list=c("mu.var", "mu.mean"))

      ### Sample tau2.mu
      tau2.shape <- tau2.gamma + n.locs/2
      tau2.rate <- tau2.phi + 0.5*t(mu - newS%*%alpha.mu)%*%(mu - newS%*%alpha.mu)
      tau2.mu <- 1/rgamma(1, shape=tau2.shape, rate=tau2.rate)

      ### Sample alpha.mu
      alpha.var <- solve((1/tau2.mu)*t(newS)%*%newS + A.prec)
      alpha.mean <- alpha.var%*%((1/tau2.mu)*t(newS)%*%mu)
      alpha.mu <- c(MASS::mvrnorm(1,alpha.mean, alpha.var))
      rm(list=c("tau2.shape", "tau2.rate", "alpha.var", "alpha.mean"))

      if ((i-burn)%%thin == 0 & i > burn) {
        mu.save[,(i-burn)/thin] <- mu
        alpha.mu.save[,(i-burn)/thin] <- alpha.mu
        tau2.mu.save[,(i-burn)/thin] <- tau2.mu
      }


    }


    ### Sample values of beta
    if (linear) {
      start = Sys.time()
      temp <- y - Jfullmu.long - Bfullxi.long - FLambda.long
      beta.var <- solve((1/sig2)*tTT + (1/tau2.beta)*Matrix::Diagonal(n=n.locs)) #ItTT + (1/tau2.beta)*Matrix::Diagonal(n=n.locs))
      beta.mean <- beta.var%*%((1/sig2)*tT%*%temp[notmissind] + (1/tau2.beta)*newS%*%alpha.beta)  #ItT%*%temp + (1/tau2.beta)*newS%*%alpha.beta)
      beta <- my_mvrnorm(beta.mean, beta.var)
      Tfullbeta.long <- Tfull%*%beta
      rm(list=c("beta.var", "beta.mean"))
      
      if(marginalize){
        precision <- kronecker(Matrix::Diagonal(n=n.locs), solve(tau2.beta*Tsub%*%t(Tsub) + sig2*diag(1, n.times)))  #Tfull <- kronecker(Matrix::Diagonal(n=n.locs), Tsub) #solve(tau2.beta*TtT + Matrix::Diagonal(x=sig2, n=dim(TtT)[1]))
        alpha.var <- solve(t(newS)%*%Matrix::t(Tfull)%*%precision%*%Tfull%*%newS + A.prec)
        alpha.mean <- alpha.var%*%t(newS)%*%Matrix::t(Tfull)%*%precision%*%temp
        alpha.beta <- as.numeric(MASS::mvrnorm(1, alpha.mean, alpha.var))
        #Tfullbeta.long <- Tfull%*%newS%*%alpha.beta
      }else{
      ### Sample alpha.beta
        alpha.var <- solve((1/tau2.beta)*t(newS)%*%newS + A.prec)
        alpha.var <- .5*alpha.var + .5*t(alpha.var)
        alpha.mean <- alpha.var%*%((1/tau2.beta)*t(newS)%*%beta)
        alpha.beta <- c(MASS::mvrnorm(1, alpha.mean, alpha.var))
      }
      
      ### Sample tau2.beta
      tau2.shape <- tau2.gamma + n.locs/2
      tau2.rate <- tau2.phi + 0.5*t(beta - newS%*%alpha.beta)%*%(beta - newS%*%alpha.beta)
      tau2.beta <- 1/rgamma(1, shape=tau2.shape, rate=tau2.rate) #scale of IG corresponds to rate of Gamma

      rm(list=c("tau2.shape", "tau2.rate", "alpha.var", "alpha.mean"))
      
      end = Sys.time()
      time.data[i,1] = end-start

      ### Save beta values
      if ((i-burn)%%thin == 0 & i > burn) {
        beta.save[,(i-burn)/thin] <- beta
        alpha.beta.save[,(i-burn)/thin] <- alpha.beta
        tau2.beta.save[,(i-burn)/thin] <- tau2.beta
      }


    }

    ### Sample Xi
    if (seasonal) {
      start = Sys.time()
      
      Tau2.xi <- Matrix::Diagonal(x=rep(1/tau2.xi, n.locs)) #Inverse of the full matrix of tau2.xi
      temp <- y - Jfullmu.long - Tfullbeta.long - FLambda.long
      
      ### Sample alpha.xi
      if(marginalize){
        #bs.basis <- mgcv::cSplineDes(doy, knots)
        #Bfull <- kronecker(Matrix::Diagonal(n=n.locs), bs.basis)
        #newS.xi <- methods::as(base::kronecker(newS, diag(n.seasn.knots)), "sparseMatrix")
        precision <- kronecker(Matrix::Diagonal(n=n.locs), solve(bs.basis%*%diag(tau2.xi)%*%t(bs.basis) + diag(sig2, nrow(bs.basis))))
        #precision <- solve(Bfull%*%diag(rep(tau2.xi, n.locs))%*%Matrix::t(Bfull) + Matrix::Diagonal(x=sig2, n=n.locs*n.times)) #n=dim(BtB)[1]))
        
        alpha.var <- solve((Matrix::t(newS.xi))%*%Matrix::t(Bfull)%*%precision%*%Bfull%*%newS.xi + methods::as(kronecker(diag(n.seasn.knots), A.prec), "sparseMatrix"))
        alpha.mean <- alpha.var%*%(Matrix::t(newS.xi))%*%Matrix::t(Bfull)%*%precision%*%temp
        alpha.xi <- as.numeric(MASS::mvrnorm(1, alpha.mean, alpha.var))
        #Bfullxi.long <- Bfullmiss%*%newS.xi%*%alpha.xi
      }else{
        alpha.var <- solve( Matrix::t(newS.xi)%*%Tau2.xi%*%newS.xi + kronecker(Matrix::Diagonal(x=1,n=n.seasn.knots), A.prec) ) #Matrix::Diagonal(x=rep(1/tau2.xi, n.locs))%*%StSI + kronecker(A.prec, Matrix::Diagonal(x=1,n=n.seasn.knots))) #Matrix::Diagonal(x=alpha.prec, n=dim(newS.xi)[2]))
        alpha.mean <- alpha.var%*%(Matrix::t(newS.xi)%*%Tau2.xi%*%xi)
        alpha.xi <- as.vector(MASS::mvrnorm(1,alpha.mean,alpha.var))
      }
      
      #Sample xi
      xi.var <- solve((1/sig2)*crossprod(Bfull) + Tau2.xi) #ItBB + (1/tau2.xi)*Matrix::Diagonal(n=n.locs*n.seasn.knots))
      xi.mean <- xi.var%*%((1/sig2)*Matrix::t(Bfull)%*%temp + Tau2.xi%*%newS.xi%*%alpha.xi) #ItB%*%temp + (1/tau2.xi)*newS.xi%*%alpha.xi)
      xi <- my_mvrnorm(xi.mean,xi.var)
      Bfullxi.long <- Bfull%*%xi
      rm(list=c("xi.var", "xi.mean"))
      
      ### Sample tau2.xi
      for(myind in 1:n.seasn.knots){
          thisxi <- seq(myind,n.seasn.knots*n.locs, by=n.seasn.knots)  
          tau2.shape <- tau2.gamma + length(thisxi)/2
          tau2.rate <- tau2.phi + 0.5 * crossprod(xi[thisxi] - newS.xi[thisxi,]%*%alpha.xi)
          tau2.xi[myind] <- 1/rgamma(1, shape=tau2.shape, rate=as.numeric(tau2.rate)) #scale of IG corresponds to rate of Gamma
      }
      
      rm(list=c("tau2.shape", "tau2.rate", "alpha.var", "alpha.mean"))
      end = Sys.time()
      time.data[i,2] = end-start

      ### Save values of xi
      if ((i-burn)%%thin == 0 & i > burn) {
        xi.save[,(i-burn)/thin] <- xi
        alpha.xi.save[,(i-burn)/thin] <- alpha.xi
        tau2.xi.save[,(i-burn)/thin] <- tau2.xi
      }


    }

    ### Sample Factors
    if (factors & i > delayFA) {
      start = Sys.time()

      ### Sample values of alphaT
      temp = y - Jfullmu.long - Tfullbeta.long - Bfullxi.long
      
      if(mean(missing)>.4){
        LamPQTmiss <- kronecker(t(Lambda.tilde), t(PQT))[,notmissind,drop=FALSE]
        lamPQTtlamPQT <- tcrossprod(LamPQTmiss) #LamPQTmiss%*%t(LamPQTmiss) 
        alphaT.var <- solve((1/sig2)*lamPQTtlamPQT + Matrix::Diagonal(x=alpha.prec, n=n.factors*n.temp.bases))
        alphaT.mean <- (1/sig2)*alphaT.var%*%LamPQTmiss%*%temp[notmissind] #matrixcalc::vec(t(PQT)%*%tempmat%*%Lambda.tilde) # matrixcalc::vec(t(QT)%*%ymat%*%Lambda.tilde) is a shortcut for t(lamQT)%*%y
      }else{
        lamPQTtlamPQT <- kronecker(t(Lambda.tilde)%*%Lambda.tilde, PQTtPQT) # This is much faster than t(kronecker(Lambda.tilde,PQT))%*%kronecker(Lambda.tilde,PQT)
        alphaT.var <- solve((1/sig2)*lamPQTtlamPQT + Matrix::Diagonal(x=alpha.prec, n=n.factors*n.temp.bases))
        tempmat = matrix(temp, nrow=n.times, ncol=n.locs)
        alphaT.mean <- (1/sig2)*alphaT.var%*%matrixcalc::vec(t(PQT)%*%tempmat%*%Lambda.tilde) # matrixcalc::vec(t(QT)%*%ymat%*%Lambda.tilde) is a shortcut for t(lamQT)%*%y
      }
      
      # alphaT <- as.vector(MASS::mvrnorm(1, alphaT.mean, alphaT.var))
      alphaT <- my_mvrnorm(alphaT.mean, alphaT.var)
      rm(list=c("alphaT.var", "alphaT.mean"))
      #Fmat = QT%*%matrix(alphaT, nrow=n.temp.bases, ncol=n.factors, byrow=FALSE)
      F.tilde = PQT%*%matrix(alphaT, nrow=n.temp.bases, ncol=n.factors, byrow=FALSE)
      end = Sys.time()
      time.data[i,3] = end-start

      ### Sample values of lambda
      start = Sys.time()

      temp = y - Jfullmu.long - Tfullbeta.long - Bfullxi.long
      #tempmat = matrix(temp, nrow=n.times, ncol=n.locs)
      PFmiss <- kronecker(Matrix::Diagonal(n=n.locs), F.tilde)[notmissind, ,drop=FALSE]
      tPFPF <- Matrix::t(PFmiss)%*%PFmiss
      #IkPFtPF <- methods::as(kronecker(diag(1, n.locs), t(F.tilde)%*%F.tilde), "sparseMatrix")
      lam.var <- solve((1/sig2)*tPFPF + Matrix::Diagonal(x=1/tau2.lambda, n=n.locs*n.factors)) #IkPFtPF + Matrix::Diagonal(x=1/tau2.lambda, n=n.locs*n.factors))
      lam.mean <- lam.var%*%((1/sig2)*(Matrix::t(PFmiss))%*%temp[notmissind] + (1/tau2.lambda)*QsI%*%alphaS) #matrixcalc::vec(t(F.tilde)%*%tempmat) + (1/tau2.lambda)*QsI%*%alphaS)

      ### which indices are fixed in the long lambda?
      Lam.index <- matrix(1:(n.factors*n.locs), nrow=n.locs, ncol=n.factors, byrow=T)
      fix.l <- c(t(Lam.index)[,factors.fixed])
      Lambda.tilde.long <- c(t(Lambda.tilde))

      lmu <- lam.mean[-fix.l] + lam.var[-fix.l, fix.l]%*%solve(lam.var[fix.l,fix.l])%*%(Lambda.tilde.long[fix.l] - lam.mean[fix.l])
      lvar <- lam.var[-fix.l, -fix.l] - lam.var[-fix.l, fix.l]%*%solve(lam.var[fix.l,fix.l])%*%lam.var[fix.l,-fix.l]
      Lambda.tilde.long[-fix.l] <- my_mvrnorm(lmu, lvar)

      Lambda.tilde <- matrix(Lambda.tilde.long, nrow=n.locs, ncol=n.factors, byrow=T)

      ### Sample values of alpha.S
      if(marginalize){
        #QsI <- methods::as(kronecker(QS, diag(1, n.factors)), "sparseMatrix")
        #QsI <- kronecker(QS, diag(1,n.factors))
        #QstQsI <- methods::as(kronecker(t(QS)%*%QS, diag(1, n.factors)), "sparseMatrix")
        IFmat <- methods::as(kronecker(diag(1,n.locs), F.tilde), "sparseMatrix")
        precision <- kronecker(Matrix::Diagonal(x=1,n=n.locs), solve(diag(sig2, n.times)+ tau2.lambda*F.tilde%*%t(F.tilde))) #solve(tau2.lambda*IFmat%*%Matrix::t(IFmat) + Matrix::Diagonal(x=sig2, n=n.times*n.locs))
        alphaS.var <- solve(Matrix::t(QsI)%*%Matrix::t(IFmat)%*%precision%*%IFmat%*%QsI + kronecker(A.lambda.prec, Matrix::Diagonal(x=1, n=n.factors)))
        alphaS.mean <- as.numeric(alphaS.var%*%Matrix::t(QsI)%*%Matrix::t(IFmat)%*%precision%*%temp)
        alphaS <- my_mvrnorm(alphaS.mean, alphaS.var)
      }else{
        alphaS.var <- solve((1/tau2.lambda)*QstQsI + kronecker(A.lambda.prec, Matrix::Diagonal(x=1, n=n.factors))) #Matrix::Diagonal(x=alpha.prec, n=n.factors*n.load.bases))
        alphaS.mean <- alphaS.var%*%((1/tau2.lambda)*matrixcalc::vec(t(Lambda.tilde)%*%QS))
        alphaS <- my_mvrnorm(alphaS.mean, alphaS.var)
      }
      rm(list=c("precision", "alphaS.var", "alphaS.mean"))

      ### Sample tau2.lambda
      tau2.shape = tau2.gamma + length(Lambda.tilde)/2
      tau2.rate = tau2.phi + 0.5*Matrix::t(Lambda.tilde.long - QsI%*%alphaS)%*%(Lambda.tilde.long - QsI%*%alphaS)
      tau2.lambda = 1/rgamma(1,shape=tau2.shape,rate=as.vector(tau2.rate))
      rm(list=c("tau2.shape", "tau2.rate"))
      end = Sys.time()
      time.data[i,4] = end-start

      FLambda.long = c(F.tilde%*%t(Lambda.tilde))

      ### Save values of FA
      if ((i-burn)%%thin == 0 & i > burn) {
        alphaT.save[,(i-burn)/thin] <- alphaT
        F.tilde.save[,(i-burn)/thin] <- matrixcalc::vec(F.tilde)
        Lambda.tilde.save[,(i-burn)/thin] <- Lambda.tilde.long
        alphaS.save[,(i-burn)/thin] <- alphaS
        tau2.lambda.save[,(i-burn)/thin] <- tau2.lambda
      }

 
    }

    ### Sample sigma2
    start=Sys.time()
    temp = y - Jfullmu.long - Tfullbeta.long - Bfullxi.long - FLambda.long
    sig2.shape = sig2.gamma + length(y[notmissind])/2
    sig2.rate = sig2.phi + 0.5*Matrix::t(as.matrix(temp[notmissind]))%*%temp[notmissind]
    sig2 = 1/rgamma(1, shape=sig2.shape, rate=as.vector(sig2.rate))
    rm(list=c("sig2.shape", "sig2.rate"))
    end=Sys.time()
    time.data[i,5] = end-start

    ### Save values of sig2
    if ((i-burn)%%thin == 0 & i > burn) {
      sig2.save[,(i-burn)/thin] <- sig2
    }



    ### Fill in missing data
    y[missind] = Jfullmu.long[missind] + Tfullbeta.long[missind] +
      Bfullxi.long[missind] + FLambda.long[missind] + rnorm(sum(missing), 0, sqrt(sig2))

    if(save.missing==T){
      if((i-burn)%%thin == 0 & i > burn){
        y.save[,(i-burn)/thin] <- y[missind]
      }
    }


    if (i %% floor(iters*.1) == 0 & verbose) {
      cat(paste("Finished iteration ", i, ": taken ", round((proc.time()[3]-start.time[3])/60,2), " minutes. \n", sep=""))
    }
    if (i == delayFA & verbose) {
      cat("Starting FA now. \n")
    }
    if (i == burn & verbose) {
      cat("Burn complete. Saving iterations now. \n")
    }
  }

  time.data$full_iter = apply(time.data,1,sum)

  if(save.time==FALSE){
    time.data <- NULL
  }

  if (verbose) cat('Finished MCMC Sampling. \n')

  output = list("mu" = coda::as.mcmc(t(mu.save)),
                "alpha.mu" = coda::as.mcmc(t(alpha.mu.save)),
                "tau2.mu" = coda::as.mcmc(t(tau2.mu.save)),
                "beta" = coda::as.mcmc(t(beta.save)),
                "alpha.beta" = coda::as.mcmc(t(alpha.beta.save)),
                "tau2.beta" = coda::as.mcmc(t(tau2.beta.save)),
                "xi" = coda::as.mcmc(t(xi.save)),
                "alpha.xi" = coda::as.mcmc(t(alpha.xi.save)),
                "tau2.xi" = coda::as.mcmc(t(tau2.xi.save)),
                "F.tilde" = coda::as.mcmc(t(F.tilde.save)),
                "alphaT" = coda::as.mcmc(t(alphaT.save)),
                "Lambda.tilde" = coda::as.mcmc(t(Lambda.tilde.save)),
                "alphaS" = coda::as.mcmc(t(alphaS.save)),
                "tau2.lambda" = coda::as.mcmc(t(tau2.lambda.save)),
                "sig2" = coda::as.mcmc(t(sig2.save)),
                "y.missing" = y.save,
                "time.data" = time.data,
                "setup.time" = setup.time,
                "model.matrices" = model.matrices,
                "factors.fixed" = factors.fixed,
                "iters" = iters,
                "y" = y,
                "ymat" = ymat,
                "missing" = missing,
                "coords" = coords,
                "doy" = doy,
                "dates" = dates,
                "knots.spatial" = knots.vec.save.spatial,
                "knots.load" = knots.vec.save.load,
                "knot.levels" = knot.levels,
                "load.style" = load.style,
                "spatial.style" = spatial.style,
                "freq.lon" = freq.lon,
                "freq.lat" = freq.lat,
                "n.times" = n.times,
                "n.locs" = n.locs,
                "n.factors" = n.factors,
                "n.seasn.knots" = n.seasn.knots,
                "n.spatial.bases" = n.spatial.bases,
                "n.temp.bases" = n.temp.bases,
                "n.load.bases" = n.load.bases,
                "draws" = dim(coda::as.mcmc(t(beta.save)))[1], 
                "mean" = mean, 
                "linear" = linear, 
                "seasonal" = seasonal, 
                "factors" = factors)

  output

}







