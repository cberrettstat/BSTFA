### Functions used to check diagnostics and convergence

#' Plot trace plots
#' @param out Output from BSTFA or BSTFAfull.
#' @param parameter Parameter to plot.  See BSTFA and BSTFAfull for parameter names.
#' @param param.range Indices of the named parameter to plot.  Default is to plot all relevant parameters.
#' @param par.mfrow Vector of length 2 indicating the number of rows and columns to divide the plotting window.
#' @param density Logical scalar indicating whether to include the density plot of the posterior draws. Default is \code{TRUE}.
#' @returns A plot containing the trace plot (and density plot when \code{density=TRUE}) of the listed parameters.
#' @author Adam Simpson
#' @examples
#' data(out.sim)
#' attach(out.sim)
#' plot_trace(out.sim, parameter='beta', param.range=1)
#' @export plot_trace
plot_trace = function(out, parameter, param.range=NULL,
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
#' @author Adam Simpson
#' @examples
#' data(utahDataList)
#' attach(utahDataList)
#' low.miss <- which(apply(is.na(TemperatureVals), 2, mean)<.02)
#' out <- BSTFA(ymat=TemperatureVals[1:50,low.miss],
#'   dates=Dates[1:50],
#'   coords=Coords[low.miss,],
#'   n.factors=2,
#'   iters=10,
#'   save.time=TRUE)
#' compute_summary(out)
#' @export compute_summary
compute_summary = function(out) {
  no.fa.iters = which(out$time.data[,3] == 0)
  fa.iters = which(out$time.data[,3] != 0)
  cat(paste('Setup Time:', round(out$setup.time,3), 'seconds. \n'))
  cat('PARAMETERS \n')
  cat(paste('Beta:', round(mean(out$time.data[,1]),3), 'seconds per iter. \n'))
  cat(paste('Xi:', round(mean(out$time.data[,2]),3), 'seconds per iter. \n'))
  cat(paste('F:', round(mean(out$time.data[fa.iters,3]),3), 'seconds per iter. \n'))
  cat(paste('Lambda:', round(mean(out$time.data[fa.iters,4]),3), 'seconds per iter. \n'))
  cat(paste('Sigma2:', round(mean(out$time.data[,5]),5), 'seconds per iter. \n'))
  cat('OVERALL PER ITERATION \n')
  cat(paste('Pre-Factor Analysis:', round(mean(out$time.data[no.fa.iters,6]),3), 'seconds. \n'))
  cat(paste('Post-Factor Analysis:', round(mean(out$time.data[fa.iters,6]),3), 'seconds. \n'))
  cat('TOTAL TIME \n')
  cat(paste('Total time: ', round(out$setup.time + apply(out$time.data,2,sum, na.rm=T)[6],3), ' seconds (', round((out$setup.time + apply(out$time.data,2,sum, na.rm=T)[6])/60,3),' minutes) for ', out$iters, ' iterations. \n', sep=""))
}


#' Check effective sample size and geweke diagnostic
#' @param out Output from BSTFA or BSTFAfull.
#' @param type Character specifying which diagnostic to compute.  Options are \code{ess} and \code{geweke}.
#' @param onlyfail Logical scalar indicating whether to save only those parameters that fail to meet the designated \code{cutoff} value. Default \code{TRUE}.
#' @param cutoff Numeric scalar indicating the cutoff value to flag parameters that haven't converged. Used when \code{onlyfail=TRUE}, ignored otherwise.
#' @returns A list containing the parameters not meeting the convergence cutoff criteria.
#' @author Adam Simpson and Candace Berrett
#' @examples
#' data(out.sim)
#' attach(out.sim)
#' convergence_diag(out.sim)
#' @importFrom coda effectiveSize
#' @importFrom coda geweke.diag
#' @export convergence_diag
convergence_diag = function(out, type='eSS', onlyfail=TRUE, cutoff=ifelse(type=='eSS',100,1.96)) {
  
  mcmcVals = list()
  
  if(onlyfail){

    if (type=='eSS') {
      if (!is.na(out$mu[1])) mcmcVals$mu = effectiveSize(out$mu)[which(effectiveSize(out$mu)<cutoff)]
      if (!is.na(out$alpha.mu[1])) mcmcVals$alpha.mu = effectiveSize(out$alpha.mu)[which(effectiveSize(out$alpha.mu)<cutoff)]
      if (!is.na(out$tau2.mu[1])) mcmcVals$tau2.mu = effectiveSize(out$tau2.mu)[which(effectiveSize(out$tau2.mu)<cutoff)]
      if (!is.na(out$beta[1])) mcmcVals$beta = effectiveSize(out$beta)[which(effectiveSize(out$beta)<cutoff)]
      if (!is.na(out$alpha.beta[1])) mcmcVals$alpha.beta = effectiveSize(out$alpha.beta)[which(effectiveSize(out$alpha.beta)<cutoff)]
      if (!is.na(out$tau2.beta[1])) mcmcVals$tau2.beta = effectiveSize(out$tau2.beta)[which(effectiveSize(out$tau2.beta)<cutoff)]
      if (!is.na(out$xi[1])) mcmcVals$xi = effectiveSize(out$xi)[which(effectiveSize(out$xi)<cutoff)]
      if (!is.na(out$alpha.xi[1])) mcmcVals$alpha.xi = effectiveSize(out$alpha.xi)[which(effectiveSize(out$alpha.xi)<cutoff)]
      if (!is.na(out$tau2.xi[1])) mcmcVals$tau2.xi = effectiveSize(out$tau2.xi)[which(effectiveSize(out$tau2.xi)<cutoff)]
      if (!is.na(out$alphaT[1])) mcmcVals$alphaT = effectiveSize(out$alphaT)[which(effectiveSize(out$alphaT)<cutoff)]
      if (!is.na(out$F.tilde[1])) mcmcVals$F.tilde = effectiveSize(out$F.tilde)[which(effectiveSize(out$F.tilde)<cutoff)]
      if (!is.na(out$Lambda[1])) mcmcVals$Lambda = effectiveSize(out$Lambda)[which(effectiveSize(out$Lambda)<cutoff)]
      if (!is.na(out$alphaS[1])) mcmcVals$alphaS = effectiveSize(out$alphaS)[which(effectiveSize(out$alphaS)<cutoff)]
      if (!is.na(out$tau2.lambda[1])) mcmcVals$tau2.lambda = effectiveSize(out$tau2.lambda)[which(effectiveSize(out$tau2.lambda)<cutoff)]
      if (!is.na(out$sig2[1])) mcmcVals$sig2 = effectiveSize(out$sig2)[which(effectiveSize(out$sig2)<cutoff)]
    }

    if (type=='geweke') {
      if (!is.na(out$mu[1])) mcmcVals$mu = geweke.diag(out$mu)[[1]][which(geweke.diag(out$mu)[[1]]>cutoff | geweke.diag(out$mu)[[1]]< -cutoff)]
      if (!is.na(out$alpha.mu[1])) mcmcVals$alpha.mu = geweke.diag(out$alpha.mu)[[1]][which(geweke.diag(out$alpha.mu)[[1]]>cutoff | geweke.diag(out$alpha.mu)[[1]]< -cutoff)]
      if (!is.na(out$tau2.mu[1])) mcmcVals$tau2.mu = geweke.diag(out$tau2.mu)[[1]][which(geweke.diag(out$tau2.mu)[[1]]>cutoff | geweke.diag(out$tau2.mu)[[1]]< -cutoff)]
      if (!is.na(out$beta[1])) mcmcVals$beta = geweke.diag(out$beta)[[1]][which(geweke.diag(out$beta)[[1]]>cutoff | geweke.diag(out$beta)[[1]]< -cutoff)]
      if (!is.na(out$alpha.beta[1])) mcmcVals$alpha.beta = geweke.diag(out$alpha.beta)[[1]][which(geweke.diag(out$alpha.beta)[[1]]>cutoff | geweke.diag(out$alpha.beta)[[1]]< -cutoff)]
      if (!is.na(out$tau2.beta[1])) mcmcVals$tau2.beta = geweke.diag(out$tau2.beta)[[1]][which(geweke.diag(out$tau2.beta)[[1]]>cutoff | geweke.diag(out$tau2.beta)[[1]]< -cutoff)]
      if (!is.na(out$xi[1])) mcmcVals$xi = geweke.diag(out$xi)[[1]][which(geweke.diag(out$xi)[[1]]>cutoff | geweke.diag(out$xi)[[1]]< -cutoff)]
      if (!is.na(out$alpha.xi[1])) mcmcVals$alpha.xi = geweke.diag(out$alpha.xi)[[1]][which(geweke.diag(out$alpha.xi)[[1]]>cutoff | geweke.diag(out$alpha.xi)[[1]]< -cutoff)]
      if (!is.na(out$tau2.xi[1])) mcmcVals$tau2.xi = geweke.diag(out$tau2.xi)[[1]][which(geweke.diag(out$tau2.xi)[[1]]>cutoff | geweke.diag(out$tau2.xi)[[1]]< -cutoff)]
      if (!is.na(out$alphaT[1])) mcmcVals$alphaT = geweke.diag(out$alphaT)[[1]][which(geweke.diag(out$alphaT)[[1]]>cutoff | geweke.diag(out$alphaT)[[1]]< -cutoff)]
      if (!is.na(out$F.tilde[1])) mcmcVals$PFmat = geweke.diag(out$F.tilde)[[1]][which(geweke.diag(out$F.tilde)[[1]]>cutoff | geweke.diag(out$F.tilde)[[1]]< -cutoff)]
      if (!is.na(out$Lambda[1])) mcmcVals$Lambda = geweke.diag(out$Lambda)[[1]][which(geweke.diag(out$Lambda)[[1]]>cutoff | geweke.diag(out$Lambda)[[1]]< -cutoff)]
      if (!is.na(out$alphaS[1])) mcmcVals$alphaS = geweke.diag(out$alphaS)[[1]][which(geweke.diag(out$alphaS)[[1]]>cutoff | geweke.diag(out$alphaS)[[1]]< -cutoff)]
      if (!is.na(out$tau2.lambda[1])) mcmcVals$tau2.lambda = geweke.diag(out$tau2.lambda)[[1]][which(geweke.diag(out$tau2.lambda)[[1]]>cutoff | geweke.diag(out$tau2.lambda)[[1]]< -cutoff)]
      if (!is.na(out$sig2[1])) mcmcVals$sig2 = geweke.diag(out$sig2)[[1]][which(geweke.diag(out$sig2)[[1]]>cutoff | geweke.diag(out$sig2)[[1]]< -cutoff)]
    }
    
    if(length(mcmcVals)==0){
      mcmcVals <- "All parameters have met the cutoff criteria."
    }
  }else{
    
    if (type=='eSS') {
      if (!is.na(out$mu[1])) mcmcVals$mu = effectiveSize(out$mu)
      if (!is.na(out$alpha.mu[1])) mcmcVals$alpha.mu = effectiveSize(out$alpha.mu)
      if (!is.na(out$tau2.mu[1])) mcmcVals$tau2.mu = effectiveSize(out$tau2.mu)
      if (!is.na(out$beta[1])) mcmcVals$beta = effectiveSize(out$beta)
      if (!is.na(out$alpha.beta[1])) mcmcVals$alpha.beta = effectiveSize(out$alpha.beta)
      if (!is.na(out$tau2.beta[1])) mcmcVals$tau2.beta = effectiveSize(out$tau2.beta)
      if (!is.na(out$xi[1])) mcmcVals$xi = effectiveSize(out$xi)
      if (!is.na(out$alpha.xi[1])) mcmcVals$alpha.xi = effectiveSize(out$alpha.xi)
      if (!is.na(out$tau2.xi[1])) mcmcVals$tau2.xi = effectiveSize(out$tau2.xi)
      if (!is.na(out$alphaT[1])) mcmcVals$alphaT = effectiveSize(out$alphaT)
      if (!is.na(out$F.tilde[1])) mcmcVals$F.tilde = effectiveSize(out$F.tilde)
      if (!is.na(out$Lambda[1])) mcmcVals$Lambda = effectiveSize(out$Lambda)
      if (!is.na(out$alphaS[1])) mcmcVals$alphaS = effectiveSize(out$alphaS)
      if (!is.na(out$tau2.lambda[1])) mcmcVals$tau2.lambda = effectiveSize(out$tau2.lambda)
      if (!is.na(out$sig2[1])) mcmcVals$sig2 = effectiveSize(out$sig2)
    }
    
    if (type=='geweke') {
      if (!is.na(out$mu[1])) mcmcVals$mu = geweke.diag(out$mu)[[1]]
      if (!is.na(out$alpha.mu[1])) mcmcVals$alpha.mu = geweke.diag(out$alpha.mu)[[1]]
      if (!is.na(out$tau2.mu[1])) mcmcVals$tau2.mu = geweke.diag(out$tau2.mu)[[1]]
      if (!is.na(out$beta[1])) mcmcVals$beta = geweke.diag(out$beta)[[1]]
      if (!is.na(out$alpha.beta[1])) mcmcVals$alpha.beta = geweke.diag(out$alpha.beta)[[1]]
      if (!is.na(out$tau2.beta[1])) mcmcVals$tau2.beta = geweke.diag(out$tau2.beta)[[1]]
      if (!is.na(out$xi[1])) mcmcVals$xi = geweke.diag(out$xi)[[1]]
      if (!is.na(out$alpha.xi[1])) mcmcVals$alpha.xi = geweke.diag(out$alpha.xi)[[1]]
      if (!is.na(out$tau2.xi[1])) mcmcVals$tau2.xi = geweke.diag(out$tau2.xi)[[1]]
      if (!is.na(out$alphaT[1])) mcmcVals$alphaT = geweke.diag(out$alphaT)[[1]]
      if (!is.na(out$F.tilde[1])) mcmcVals$F.tilde = geweke.diag(out$F.tilde)[[1]]
      if (!is.na(out$Lambda[1])) mcmcVals$Lambda = geweke.diag(out$Lambda)[[1]]
      if (!is.na(out$alphaS[1])) mcmcVals$alphaS = geweke.diag(out$alphaS)[[1]]
      if (!is.na(out$tau2.lambda[1])) mcmcVals$tau2.lambda = geweke.diag(out$tau2.lambda)[[1]]
      if (!is.na(out$sig2[1])) mcmcVals$sig2 = geweke.diag(out$sig2)[[1]]
    }
  }
  mcmcVals

}

#' Compute log-likelihood
#' @param out Output from BSTFA or BSTFAfull.
#' @param verbose Logical scalar indicating whether to print status of the log-likelihood computation.  Default is \code{FALSE}.
#' @param addthin Numeric scalar indicating the number of additional draws to thin by to reduce the computation time.  Default is \code{1} (no additional thinning).
#' @returns A matrix of size \code{n.times*n.locs} by \code{draws} log-likelihood values for each observation and each posterior draw.
#' @author Adam Simpson and Candace Berrett
#' @examples
#' data(out.sim)
#' attach(out.sim)
#' loglik <- computeLogLik(out.sim, addthin=2)
#' @export computeLogLik
computeLogLik <- function(out, verbose=FALSE, addthin=1) {
  y = out$y
  mu = predictBSTFA(out=out,
                    type='all', pred.int=FALSE)
  myseq <- seq(1,out$draws,by=addthin)
  log_lik = matrix(0,nrow=out$n.times*out$n.locs,
                          ncol=length(myseq))
  if (verbose) cat('Starting Log-likelihood calculation \n')
  for (d in 1:length(myseq)) {
      log_lik[,d] = dnorm(y,mu[,myseq[d]],sd=sqrt(out$sig2[myseq[d]]),log=TRUE)
    if (verbose & (d %% floor(length(myseq)*.1) == 0)){ cat(paste('Draw', d, '\n')) }
  }
  log_lik
}


