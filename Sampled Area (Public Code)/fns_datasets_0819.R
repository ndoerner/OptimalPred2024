
# This code is modified from fns_datasets_0805.R
## correct the mistake in vcov.jack function use samparea!!! Not sampnum!!!
### vcov.jack.parallel  lapply -> mclapply

##################################i###########################
####################informative sampling######################
##############################################################

#######data generation

source("fns_population_gen_sampling_0623.R")

sseexppi <- function(coef, ypi, y, x, areaspi){
  
  kipi <-  coef[ 2 + (1:length(unique(areaspi)))]
  
  # areaspi is a numeric vector. names(kipi) is a character vector of it.
  names(kipi) <- sort(unique(areaspi))
  sum( (ypi - kipi[as.character(areaspi)]*exp( coef[1]*x+coef[2]*y ))^2) # VVV
  
}

isSingular2 <-function (x, tol = 1e-05) 
{
  lwr <- getME(x, "lower")
  theta <- getME(x, "theta")
  any(theta[lwr == 0] < tol)
}


P.LikHood = function(s2u,datsamp,lmersamp){
  
  #  
  sig2ehatsamp <- summary(lmersamp)$sigma^2
  samparea = unique(datsamp$area)
  nisneword = tapply(datsamp$area, datsamp$area, length )
  
  #Sigma_i
  block = lapply(nisneword, function(ni){
    Sigmai = s2u*matrix(rep(1,ni^2),nrow=ni,ncol=ni)+diag(rep(sig2ehatsamp,ni))
  })
  
  #det(Sigma_i) & log(Sigam) in L_p(A)
  Sigma= sapply(1:length(nisneword),function(i){det(block[[i]])})
  SigDet = sum(log(Sigma))
  
  #(Simga_i)^{-1}
  SigInv = lapply(1:length(nisneword),function(i){solve(block[[i]])})
  
  #P.comp = X'Sigma_{i}^{-1}X and yPiy 
  
  X= cbind(rep(1,nrow(datsamp)),datsamp$x)
  LargeSigma = bdiag(block);LargeSigma = as.matrix(LargeSigma)
  ILargeSigma = solve(LargeSigma)
  P.comp = t(X)%*%ILargeSigma%*%X
  P = ILargeSigma - ILargeSigma%*%X%*%solve(P.comp)%*%t(X)%*%ILargeSigma 
  
  
  logLikAdj = log(s2u)-0.5*SigDet-t(datsamp$yij)%*%P%*%datsamp$yij/2
  return(-logLikAdj)
  
}

## Estimate the parameters by method
gen.est.par = function(datsamp, par_method = "Ori"){
  
  ## PS General
  lmersamp <- lmer(yij~x + (1|area), data = datsamp,REML=FALSE ) #f(yij|i=1,j=1)
  if(isSingular2(lmersamp)){
    #Adjusted Likelihood
    sig2ehatsamp <- summary(lmersamp)$sigma^2
    sig2uhatsamp <- optim(par=VarCorr(lmersamp)[[1]][1],fn=P.LikHood,datsamp=datsamp,method="Brent",lower=0,upper=200,lmersamp=lmersamp)$par
  }else{
    ##### Implement EBUP simulation method:
    sig2ehatsamp <- summary(lmersamp)$sigma^2
    sig2uhatsamp <- VarCorr(lmersamp)[[1]][1]
  }
  
  ################################### CHANGED for JACKKNIFE ##########################################################################
  lmpireg <- lm(log(1/datsamp$probij )~ datsamp$x + datsamp$yij + as.factor(datsamp$area)-1)
  areaspi <- datsamp$area
  lmpireg$coef[ 2 + (1:length(unique(areaspi)))] <- exp(lmpireg$coef[ 2 + (1:length(unique(areaspi)))])
  
  o=tryCatch(expr =    nls(1/probij~kappa[as.factor(area)]*exp(b*yij+a*x),
                           data=datsamp,
                           start =list(a=lmpireg$coef[1],b=lmpireg$coef[2],kappa=lmpireg$coef[2+1:length(unique(areaspi))])),
             error = function(e){TRUE}
  )
  
  if(is.list(o)==TRUE){
    coefunit <- optim(summary(o)$parameters[,1], sseexppi, ypi = 1/datsamp$probij ,y = datsamp$y, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 100000))$par
  } else {
    coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 10000))$par
  }
  
  # coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 10000))$par
  ######################################################################################################################################################################
  
  
  
  # Pseudo EBP function
  samp.wij = 1/datsamp$probij
  samp.norm.wij =  (1/datsamp$probij)/tapply(1/datsamp$probij,datsamp$area,sum)[as.character(datsamp$area)]
  
  ybaris.w = tapply(samp.norm.wij*datsamp$yij, datsamp$area, sum)
  xbaris.w = tapply(samp.norm.wij*datsamp$x, datsamp$area,sum)
  deltai2 = tapply(samp.norm.wij^2, datsamp$area,sum)
  
  gammahatis.w = sig2uhatsamp/(sig2uhatsamp + sig2ehatsamp*deltai2)
  
  
  if(par_method == "Ori"){
    betahatsamp <- fixef(lmersamp)
    nis <- tapply(datsamp$area, datsamp$area, length ) # should be length!!!!!!!!!!! not sum!!!!!!!!!!!
    ybaris <- tapply(datsamp$yij, datsamp$area, mean)
    xbaris <- tapply(datsamp$x, datsamp$area, mean)
    gammahatis <- sig2uhatsamp/(sig2uhatsamp + sig2ehatsamp/nis)
    uhatis <- gammahatis*(ybaris - betahatsamp[1] - betahatsamp[2]*xbaris)
    vcondus <- gammahatis*sig2ehatsamp/nis
    
    theta.hat =  list(betahat = betahatsamp
                      , sig2uhat = sig2uhatsamp
                      , sig2ehat = sig2ehatsamp
                      , xbar = xbaris
                      , ybar = ybaris
                      , gammahat = gammahatis 
                      , coefunit_12 = coefunit[1:2]
    ) 
    
  } else if (par_method == "PS"){
    
    betahatsamp <- fixef(lmersamp)
    coefunit[2] <- 0
    
    theta.hat =  list(betahat = betahatsamp
                      , sig2uhat = sig2uhatsamp
                      , sig2ehat = sig2ehatsamp
                      , xbar = xbaris.w
                      , ybar = ybaris.w
                      , gammahat = gammahatis.w 
                      , coefunit_12 = coefunit[1:2]
    ) 
    
  }  else if (par_method == "AW"){
    
    betahatsamp.aw.d = lapply(1:length(xbaris.w), function(.){gammahatis.w[.]*(c(1,xbaris.w[.])%*%t(c(1,xbaris.w[.])))}) %>% Reduce("+",.)
    betahatsamp.aw.n = lapply(1:length(xbaris.w), function(.){gammahatis.w[.]*c(1,xbaris.w[.])*ybaris.w[.]}) %>% Reduce("+",.)
    betahatsamp.aw = solve(betahatsamp.aw.d)%*%betahatsamp.aw.n
    coefunit[2] = 0
    
    
    theta.hat =  list(betahat = betahatsamp.aw
                      , sig2uhat = sig2uhatsamp
                      , sig2ehat = sig2ehatsamp
                      , xbar = xbaris.w
                      , ybar = ybaris.w
                      , gammahat = gammahatis.w 
                      , coefunit_12 = coefunit[1:2]
    )
    
  }else if (par_method == "W") {
    x.w.new = cbind((1-gammahatis.w)[as.character(datsamp$area)], datsamp$x-(gammahatis.w*xbaris.w)[as.character(datsamp$area)])
    betahatsamp.w.d = lapply(1:nrow(datsamp), function(.){samp.wij[.]*(c(1,datsamp$x[.])%*%t(x.w.new[.,]))}) %>% Reduce("+",.)
    betahatsamp.w.n = lapply(1:nrow(datsamp), function(.){samp.wij[.]*x.w.new[.,]*datsamp$yij[.]}) %>% Reduce("+",.)
    betahatsamp.w = solve(betahatsamp.w.d)%*%betahatsamp.w.n
    coefunit[2] = 0
    
    theta.hat =  list(betahat = betahatsamp.w
                      , sig2uhat = sig2uhatsamp
                      , sig2ehat = sig2ehatsamp
                      , xbar = xbaris.w
                      , ybar = ybaris.w
                      , gammahat = gammahatis.w 
                      , coefunit_12 = coefunit[1:2]
    )
  } else {
    print("Be careful your method")
  }
  
  return(theta.hat)
  
}


gen.vcov.jack= function(datsamp, par_method){ ## has been changed
  
  sampnum = length(unique(datsamp$area))
  samparea = unique(datsamp$area)   # updated 08162021!!!!!!!
  
  if (par_method == "Ori"){
    
    par.mat = lapply(1:sampnum, function(b){
      
      res = gen.est.par(datsamp[datsamp$area!=samparea[b],], par_method)   ### changed from datsamp[datsamp$area!=b,]
      # estpar
      c(res$betahat, res$sig2uhat, res$sig2ehat, res$coefunit_12) %>% 
        setNames(c("betahat0","betahat1","sig2uhat","sig2ehat","coef_x","coef_yij")) %>% 
        return(.)
      
    }) %>% do.call("rbind",.)
    
  } else {
    
    par.mat = lapply(1:sampnum, function(b){
      
      res = gen.est.par(datsamp[datsamp$area!=samparea[b],], par_method)
      # estpar
      c(res$betahat, res$sig2uhat, res$sig2ehat) %>% 
        setNames(c("betahat0","betahat1","sig2uhat","sig2ehat")) %>%
        return(.) 
      
    }) %>% do.call("rbind",.)
    
  }
  
  par.mean = apply(par.mat,2,mean)
  par.sum = lapply(1:sampnum, function(i){
    
    diff = as.matrix( par.mat[i,]-par.mean)
    diff.mat = diff%*%t(diff)
    
  })
  
  
  return((sampnum-1)* Reduce("+",par.sum)/sampnum)
}





gen.vcov.jack.parallel = function(datsamp, par_method, core_num){ ## has been changed
  
  sampnum = length(unique(datsamp$area))
  samparea = unique(datsamp$area)   # updated 08162021!!!!!!!
  
  if (par_method == "Ori"){
    
    par.mat = mclapply(1:sampnum, function(b){
      
      res = gen.est.par(datsamp[datsamp$area!=samparea[b],], par_method)   ### changed from datsamp[datsamp$area!=b,]
      # estpar
      c(res$betahat, res$sig2uhat, res$sig2ehat, res$coefunit_12) %>% 
        setNames(c("betahat0","betahat1","sig2uhat","sig2ehat","coef_x","coef_yij")) %>% 
        return(.)

    }, mc.cores = core_num) %>% do.call("rbind",.)

  } else {
    
    par.mat = lapply(1:sampnum, function(b){
      
      res = gen.est.par(datsamp[datsamp$area!=samparea[b],], par_method)
      # estpar
      c(res$betahat, res$sig2uhat, res$sig2ehat) %>% 
        setNames(c("betahat0","betahat1","sig2uhat","sig2ehat")) %>%
        return(.) 
      
      }) %>% do.call("rbind",.)

  }
  
  par.mean = apply(par.mat,2,mean)
  par.sum = lapply(1:sampnum, function(i){
    
    diff = as.matrix( par.mat[i,]-par.mean)
    diff.mat = diff%*%t(diff)
    
  })
  
  
  return((sampnum-1)* Reduce("+",par.sum)/sampnum)
}

dat.gen = function(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, par_method = "Ori",stratindarea, stratindpop, core_num){
  
  # Finite population Generation
  ranpopPS <- genpopPS(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij)
  # ranpopPS <- genpopPSNoTrunc(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij)
  usD <-ranpopPS[[1]]; eij <- ranpopPS[[2]]; yij <- ranpopPS[[3]]
  
  sampPS <- gensampPS(D, usD, eij, yij, sig2u, sig2e, stratindarea, stratindpop, Nis,areafacpop) 
  probi <- sampPS[[1]]; Iareas <- sampPS[[2]]; samparea <- sampPS[[3]]; 
  probij <- sampPS[[4]]; Iunits <- sampPS[[5]]
  
  names(probi) <- names(Iareas) <- 1:D
  
  # Make data sets
  
  datpop <- data.frame(area = areafacpop, strat = stratindpop, x = xij, 
                       probi = probi[areafacpop], probij = probij, yij = yij, Iareas = Iareas[areafacpop], Iunits = Iunits )
  
  datsamparea <- subset(datpop, Iareas == 1)
  datsamparea$ID <- 1:nrow(datsamparea)
  
  datnonsamparea <-subset(datpop, Iareas==0)
  
  datsamp  <- subset(datpop, Iareas == 1 & Iunits == 1)
  
  ##  Fit sample model to sampled data:
  
  lmersamp <- lmer(yij~x + (1|area), data = datsamp,REML=FALSE ) #f(yij|i=1,j=1)

  sig2ehatsamp <- summary(lmersamp)$sigma^2  
  
  if(isSingular2(lmersamp)){
    #Adjusted Likelihood
    sig2uhatsamp <- optim(par=VarCorr(lmersamp)[[1]][1],fn=P.LikHood,datsamp=datsamp,method="Brent",lower=0,upper=200,lmersamp=lmersamp)$par
    
  }else{
    sig2uhatsamp <- VarCorr(lmersamp)[[1]][1]
  }
  
  # ##### Implement EBUP simulation method:

  betahatsamp <- fixef(lmersamp)
  nis <- tapply(datsamp$area, datsamp$area, length ) # should be length!!!!!!!!!!! not sum!!!!!!!!!!!
  n <-sum(nis)
  ybaris <- tapply(datsamp$yij, datsamp$area, mean)
  xbaris <- tapply(datsamp$x, datsamp$area, mean)
  gammahatis <- sig2uhatsamp/(sig2uhatsamp + sig2ehatsamp/nis)
  uhatis <- gammahatis*(ybaris - betahatsamp[1] - betahatsamp[2]*xbaris)
  vcondus <- gammahatis*sig2ehatsamp/nis
  areaspi <- datsamp$area
  
  #### estimated unit weight model parameters  ### CHANGED PART IT WORKS WELL for JACKKNIFE
  lmpireg <- lm(log(1/datsamp$probij )~ datsamp$x + datsamp$yij + as.factor(datsamp$area)-1)
  areaspi <- datsamp$area
  lmpireg$coef[ 2 + (1:length(unique(areaspi)))] <- exp(lmpireg$coef[ 2 + (1:length(unique(areaspi)))])
  
  o=tryCatch(expr =    nls(1/probij~kappa[as.factor(area)]*exp(b*yij+a*x),
                           data=datsamp,
                           start =list(a=lmpireg$coef[1],b=lmpireg$coef[2],kappa=lmpireg$coef[2+1:length(unique(areaspi))])),
             error = function(e){TRUE}
  )
  
  if(is.list(o)==TRUE){
    coefunit <- optim(summary(o)$parameters[,1], sseexppi, ypi = 1/datsamp$probij ,y = datsamp$y, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 100000))$par
  } else {
    coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 10000))$par
  }
  # coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 10000))$par
  
  #update 02092021
  # vcov.jack = gen.vcov.jack.parallel(datsamp, par_method, core_num)  ## changed 08212021
  vcov.jack = gen.vcov.jack(datsamp, par_method)  ## changed 08212021
  
  ############################################################################
  # Pseudo EBP function
  samp.wij = 1/datsamp$probij
  samp.norm.wij =  (1/datsamp$probij)/tapply(1/datsamp$probij,datsamp$area,sum)[as.character(datsamp$area)]
  
  ybaris.w = tapply(samp.norm.wij*datsamp$yij, datsamp$area, sum)
  xbaris.w = tapply(samp.norm.wij*datsamp$x, datsamp$area,sum)
  deltai2 = tapply(samp.norm.wij^2, datsamp$area,sum)
  
  gammahatis.w = sig2uhatsamp/(sig2uhatsamp + sig2ehatsamp*deltai2)
  
  
  betahatsamp.aw.d = lapply(1:length(xbaris.w), function(.){gammahatis.w[.]*(c(1,xbaris.w[.])%*%t(c(1,xbaris.w[.])))}) %>% Reduce("+",.)
  betahatsamp.aw.n = lapply(1:length(xbaris.w), function(.){gammahatis.w[.]*c(1,xbaris.w[.])*ybaris.w[.]}) %>% Reduce("+",.)
  
  betahatsamp.aw = solve(betahatsamp.aw.d)%*%betahatsamp.aw.n
  
  
  x.w.new = cbind((1-gammahatis.w)[as.character(datsamp$area)], datsamp$x-(gammahatis.w*xbaris.w)[as.character(datsamp$area)])
  
  betahatsamp.w.d = lapply(1:nrow(datsamp), function(.){
    samp.wij[.]*(c(1,datsamp$x[.])%*%t(x.w.new[.,])) }) %>% Reduce("+",.)
  
  betahatsamp.w.n = lapply(1:nrow(datsamp), function(.){
    samp.wij[.]*x.w.new[.,]*datsamp$yij[.] }) %>% Reduce("+",.)
  
  betahatsamp.w = solve(betahatsamp.w.d)%*%betahatsamp.w.n

  if(par_method == "Ori"){
    par.mu = c(betahatsamp, sig2uhatsamp, sig2ehatsamp, coefunit[1], coefunit[2])
  }else if (par_method == "PS"){
    par.mu = c(betahatsamp, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0)
  }else if (par_method == "AW"){
    par.mu = c(betahatsamp.aw, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0)
  }else if (par_method == "W"){
    par.mu = c(betahatsamp.w, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0)
  }else {par.mu = "Wrong"; print(par.mu)}
  
  
  ###################################################################
  
  
  list(datpop = datpop, datsamparea = datsamparea, datsamp = datsamp, datnonsamparea = datnonsamparea,samparea = samparea,
       Iareas = Iareas, Iunits = Iunits, sig2uhatsamp = sig2uhatsamp, sig2ehatsamp = sig2ehatsamp, betahatsamp = betahatsamp, 
       nis = nis, ybaris= ybaris, xbaris=xbaris,gammahatis = gammahatis,  uhatis = uhatis, vcondus = vcondus, areaspi = areaspi, 
       coefunit = coefunit
       , vcov.jack = vcov.jack
       , par.mu = par.mu
       , samp.norm.wij =  samp.norm.wij
       , betahatsamp.aw = betahatsamp.aw
       , betahatsamp.w = betahatsamp.w
       , deltai2 = deltai2
       , gammahatis.w = gammahatis.w
       , ybaris.w = ybaris.w
       , xbaris.w = xbaris.w
  )

}
