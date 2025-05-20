### Functions used to check diagnostics and convergence

#' Plot trace plots
#' @param out Output from BSTFA or BSTFAfull.
#' @param parameter Parameter to plot.  See BSTFA and BSTFAfull for parameter names.
#' @param param.range Indices of the named parameter to plot.  Default is to plot all relevant parameters.
#' @param par.mfrow Vector of length 2 indicating the number of rows and columns to divide the plotting window.
#' @param density Logical scalar indicating whether to include the density plot of the posterior draws. Default is \code{TRUE}.
#' @returns A plot containing the trace plot (and density plot when \code{density=TRUE}) of the listed parameters.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' plot.trace(out, parameter="beta", param.range=1)
#' @export plot.trace
plot.trace = function(out, parameter, param.range=NULL,
                      par.mfrow=c(1,1), density=TRUE) {

  vals = out[[parameter]]

  if (is.null(param.range)) ind=1:dim(vals)[2]
  else ind=param.range

  if (density==TRUE) {
    for (i in ind) {
      plot(vals[,i], type='l', main = paste(parameter, i))
    }
  } else {
    for (i in ind) {
      plot(vals[,i], type='l', main = paste(parameter, i),
           density=FALSE)
    }
  }
}

#' Print computation summary
#' @param out Output from BSTFA or BSTFAfull.
#' @returns Prints the computation time per iteration for each parameter.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' computation.summary(out)
#' @export computation.summary
computation.summary = function(out) {
  no.fa.iters = which(out$time.data[,3] == 0)
  fa.iters = which(out$time.data[,3] != 0)
  print(paste('Setup Time:', round(out$setup.time,3), 'seconds.'))
  print('PARAMETERS')
  print(paste('Beta:', round(mean(out$time.data[,1]),3), 'seconds per iter.'))
  print(paste('Xi:', round(mean(out$time.data[,2]),3), 'seconds per iter.'))
  print(paste('F:', round(mean(out$time.data[fa.iters,3]),3), 'seconds per iter.'))
  print(paste('Lambda:', round(mean(out$time.data[fa.iters,4]),3), 'seconds per iter.'))
  print(paste('Sigma2:', round(mean(out$time.data[,5]),5), 'seconds per iter.'))
  print('OVERALL PER ITERATION')
  print(paste('Pre-Factor Analysis:', round(mean(out$time.data[no.fa.iters,6]),3), 'seconds.'))
  print(paste('Post-Factor Analysis:', round(mean(out$time.data[fa.iters,6]),3), 'seconds.'))
  print('TOTAL TIME')
  print(paste('Total time: ', round(out$setup.time + apply(out$time.data,2,sum)[6],3), ' seconds (', round((out$setup.time + apply(out$time.data,2,sum)[6])/60,3),' minutes) for ', out$iters, ' iterations.', sep=""))
}


#' Check effective sample size and geweke diagnostic
#' @param out Output from BSTFA or BSTFAfull.
#' @param type Character specifying which diagnostic to compute.  Options are \code{ess} and \code{geweke}.
#' @param cutoff Numeric scalar indicating the cutoff value to flag parameters that haven't converged.
#' @returns A list containing convergence diagnostic for all parameters.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' check.convergence(out)
#' @importFrom coda effectiveSize
#' @importFrom coda geweke.diag
#' @export check.convergence
check.convergence = function(out, type='eSS', cutoff=ifelse(type=='eSS',100,0.001)) {

  mcmcVals = list()
  if (type=='eSS') {
    if (sum(out$mu)!=0) mcmcVals$mu = effectiveSize(out$mu)[which(effectiveSize(out$mu)<cutoff)]
    if (sum(out$alpha.mu)!=0) mcmcVals$alpha.mu = effectiveSize(out$alpha.mu)[which(effectiveSize(out$alpha.mu)<cutoff)]
    if (sum(out$tau2.mu)!=0) mcmcVals$tau2.mu = effectiveSize(out$tau2.mu)[which(effectiveSize(out$tau2.mu)<cutoff)]
    if (sum(out$beta)!=0) mcmcVals$beta = effectiveSize(out$beta)[which(effectiveSize(out$beta)<cutoff)]
    if (sum(out$alpha.beta)!=0) mcmcVals$alpha.beta = effectiveSize(out$alpha.beta)[which(effectiveSize(out$alpha.beta)<cutoff)]
    if (sum(out$tau2.beta)!=0) mcmcVals$tau2.beta = effectiveSize(out$tau2.beta)[which(effectiveSize(out$tau2.beta)<cutoff)]
    if (sum(out$xi)!=0) mcmcVals$xi = effectiveSize(out$xi)[which(effectiveSize(out$xi)<cutoff)]
    if (sum(out$alpha.xi)!=0) mcmcVals$alpha.xi = effectiveSize(out$alpha.xi)[which(effectiveSize(out$alpha.xi)<cutoff)]
    if (sum(out$tau2.xi)!=0) mcmcVals$tau2.xi = effectiveSize(out$tau2.xi)[which(effectiveSize(out$tau2.xi)<cutoff)]
    if (sum(out$alphaT)!=0) mcmcVals$alphaT = effectiveSize(out$alphaT)[which(effectiveSize(out$alphaT)<cutoff)]
    if (sum(out$PFmat)!=0) mcmcVals$PFmat = effectiveSize(out$PFmat)[which(effectiveSize(out$PFmat)<cutoff)]
    if (sum(out$Lambda)!=0) mcmcVals$Lambda = effectiveSize(out$Lambda)[which(effectiveSize(out$Lambda)<cutoff)]
    if (sum(out$alphaS)!=0) mcmcVals$alphaS = effectiveSize(out$alphaS)[which(effectiveSize(out$alphaS)<cutoff)]
    if (sum(out$tau2.lambda)!=0) mcmcVals$tau2.lambda = effectiveSize(out$tau2.lambda)[which(effectiveSize(out$tau2.lambda)<cutoff)]
    if (sum(out$sig2)!=0) mcmcVals$sig2 = effectiveSize(out$sig2)[which(effectiveSize(out$sig2)<cutoff)]
  }

  if (type=='geweke') {
    if (sum(out$mu)!=0) mcmcVals$mu = geweke.diag(out$mu)[[1]][which(geweke.diag(out$mu)[[1]]>cutoff | geweke.diag(out$mu)[[1]]< -cutoff)]
    if (sum(out$alpha.mu)!=0) mcmcVals$alpha.mu = geweke.diag(out$alpha.mu)[[1]][which(geweke.diag(out$alpha.mu)[[1]]>cutoff | geweke.diag(out$alpha.mu)[[1]]< -cutoff)]
    if (sum(out$tau2.mu)!=0) mcmcVals$tau2.mu = geweke.diag(out$tau2.mu)[[1]][which(geweke.diag(out$tau2.mu)[[1]]>cutoff | geweke.diag(out$tau2.mu)[[1]]< -cutoff)]
    if (sum(out$beta)!=0) mcmcVals$beta = geweke.diag(out$beta)[[1]][which(geweke.diag(out$beta)[[1]]>cutoff | geweke.diag(out$beta)[[1]]< -cutoff)]
    if (sum(out$alpha.beta)!=0) mcmcVals$alpha.beta = geweke.diag(out$alpha.beta)[[1]][which(geweke.diag(out$alpha.beta)[[1]]>cutoff | geweke.diag(out$alpha.beta)[[1]]< -cutoff)]
    if (sum(out$tau2.beta)!=0) mcmcVals$tau2.beta = geweke.diag(out$tau2.beta)[[1]][which(geweke.diag(out$tau2.beta)[[1]]>cutoff | geweke.diag(out$tau2.beta)[[1]]< -cutoff)]
    if (sum(out$xi)!=0) mcmcVals$xi = geweke.diag(out$xi)[[1]][which(geweke.diag(out$xi)[[1]]>cutoff | geweke.diag(out$xi)[[1]]< -cutoff)]
    if (sum(out$alpha.xi)!=0) mcmcVals$alpha.xi = geweke.diag(out$alpha.xi)[[1]][which(geweke.diag(out$alpha.xi)[[1]]>cutoff | geweke.diag(out$alpha.xi)[[1]]< -cutoff)]
    if (sum(out$tau2.xi)!=0) mcmcVals$tau2.xi = geweke.diag(out$tau2.xi)[[1]][which(geweke.diag(out$tau2.xi)[[1]]>cutoff | geweke.diag(out$tau2.xi)[[1]]< -cutoff)]
    if (sum(out$alphaT)!=0) mcmcVals$alphaT = geweke.diag(out$alphaT)[[1]][which(geweke.diag(out$alphaT)[[1]]>cutoff | geweke.diag(out$alphaT)[[1]]< -cutoff)]
    if (sum(out$PFmat)!=0) mcmcVals$PFmat = geweke.diag(out$PFmat)[[1]][which(geweke.diag(out$PFmat)[[1]]>cutoff | geweke.diag(out$PFmat)[[1]]< -cutoff)]
    if (sum(out$Lambda)!=0) mcmcVals$Lambda = geweke.diag(out$Lambda)[[1]][which(geweke.diag(out$Lambda)[[1]]>cutoff | geweke.diag(out$Lambda)[[1]]< -cutoff)]
    if (sum(out$alphaS)!=0) mcmcVals$alphaS = geweke.diag(out$alphaS)[[1]][which(geweke.diag(out$alphaS)[[1]]>cutoff | geweke.diag(out$alphaS)[[1]]< -cutoff)]
    if (sum(out$tau2.lambda)!=0) mcmcVals$tau2.lambda = geweke.diag(out$tau2.lambda)[[1]][which(geweke.diag(out$tau2.lambda)[[1]]>cutoff | geweke.diag(out$tau2.lambda)[[1]]< -cutoff)]
    if (sum(out$sig2)!=0) mcmcVals$sig2 = geweke.diag(out$sig2)[[1]][which(geweke.diag(out$sig2)[[1]]>cutoff | geweke.diag(out$sig2)[[1]]< -cutoff)]
  }

  mcmcVals

}

#' Compute log-likelihood
#' @param out Output from BSTFA or BSTFAfull.
#' @param verbose Logical scalar indicating whether to print status of the log-likelihood computation.  Default is \code{FALSE}.
#' @param addthin Numeric scalar indicating the number of additional draws to thin by to reduce the computation time.  Default is \code{1} (no additional thinning).
#' @returns A matrix of size \code{n.times*n.locs} by \code{draws} log-likelihood values for each observation and each posterior draw.
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals, dates=Dates, coords=Coords, iters=100)
#' loglik <- computeLogLik(out, addthin=2)
#' @export computeLogLik
computeLogLik <- function(out, verbose=FALSE, addthin=1) {
  y = out$y
  mu = predictBSTFA(out=out,
                    type='all')
  myseq <- seq(1,out$draws,by=addthin)
  log_lik = matrix(0,nrow=out$n.times*out$n.locs,
                          ncol=length(myseq))
  if (verbose) print('Starting Log-likelihood calculation')
  for (d in 1:length(myseq)) {
      log_lik[,d] = dnorm(y,mu[,myseq[d]],sd=out$sig2[myseq[d]],log=TRUE)
    if (verbose) print(paste('Draw', d))
  }
  log_lik
}




# compute.DIC <- function(out){
#   #in-sample DIC only (this function removes ALL data treated as missing)
#
#   n.total <- out$n.times*out$n.locs
#
#   if(sum(out$mu)!=0){
#     meanmat <- kronecker(diag(1, out$n.locs), rep(1, out$n.times))
#     mu.mean <- apply(out$mu, 2, mean)
#     mulong <- mean.mat%*%mu.mean
#   }else{
#     mulong <- rep(0, n.total)
#   }
#
#   if(sum(out$beta)!=0){
#     Tfull <- kronecker(diag(1, out$n.locs), out$model.matrices$linear.Tsub)
#     beta.mean <- apply(out$beta, 2, mean)
#     betalong <- Tfull%*%beta.mean
#   }else{
#     betalong <- rep(0, n.total)
#   }
#
#   if(sum(out$xi)!=0){
#     knots <- seq(1, 366, length=out$n.seasn.knots+1)
#     bs.basis <- cSplineDes(doy, knots)
#     Bfull <- kronecker(diag(1, out$n.locs), out$model.matrices$seasonal.bs.basis)
#     xi.mean <- apply(out$xi, 2, mean)
#     xilong <- Bfull%*%xi.mean
#   }else{
#     xilong <- rep(0, n.total)
#   }
#
#   if(sum(out$Lambda)!=0){
#     PF.mat.mean <- matrix(apply(out$PFmat, 2, mean),nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)
#     Lambda.mean <- matrix(apply(out$Lambda, 2, mean),nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)
#     factors.mean <- PF.mat.mean%*%t(Lambda.mean)
#     factorslong <- c(factors.mean)
#   }else{
#     factorslong <- rep(0, n.total)
#   }
#
#   sigma2.mean <- mean(out$sig2)
#
#   D.thetabar <- sum(-2*dnorm(out$y[!out$missing], (mulong+betalong+xilong+factorslong)[!out$missing], sqrt(sigma2.mean), log=T), na.rm=T)
#   D <- rep(0, draws)
#   for(ii in 1:length(D)){
#     this.mean <- 0
#     if(sum(out$mu)!=0){
#       this.mean <- this.mean + mean.mat%*%out$mu[ii,]
#     }
#     if(sum(out$beta)!=0){
#       this.mean <- this.mean + Tfull%*%out$beta[ii,]
#     }
#     if(sum(out$xi)!=0){
#       this.mean <- this.mean + Bfull%*%out$xi[ii,]
#     }
#     if(sum(out$Lambda)!=0){
#       this.mean <- this.mean + c(matrix(out$PFmat[ii,],nrow=out$n.times,ncol=out$n.factors,byrow=FALSE)%*%t(matrix(out$Lambda[ii,],nrow=out$n.locs,ncol=out$n.factors,byrow=TRUE)))
#     }
#     D[ii] <- sum(-2*dnorm(out$y[!out$missing], this.mean[!out$missing], sqrt(out$sig2), log=T), na.rm=T)
#   }
#   D.bar <- mean(D)
#   DIC <- 2*D.bar - D.thetabar
#   DIC
# }
