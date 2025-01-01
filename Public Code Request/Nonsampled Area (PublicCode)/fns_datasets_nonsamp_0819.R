
# This code is modified from fns_datasets_0623.R
## parweight added for nonsampled areas : this is for the optimal function to estimate area-level probability prob_i model parameters.
## coefarea added in dat.gen

# This code is modified from fns_datasets_nonsamp_0805.R
## correct the mistake in vcov.jack function use samparea!!! Not sampnum!!!


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
Fct <- function(coef,W,U,log_pi, uhatis, vcondus){
  
  par1 = coef[1]
  par2 = coef[2]
  tau = abs(coef[3]); itau = 1/tau
  
  samp.n = length(vcondus)
  Integral_i = sapply(1:samp.n, function(.){
    f_x = exp( -((-log_pi[.]-par1-par2*(sqrt(2)*vcondus[.]*U+uhatis[.]))/(sqrt(2)*tau))^2)*itau/sqrt(2)/pi
    log(sum(W*f_x))
    
  })
  
  
  Neg_ll = - sum(Integral_i)
  
  return(Neg_ll)
}
# Gaussian Hermite method is used.
Coef_Area_GH = function(datsamp,lmersamp, M = 100, uhatis, vcondus){
  # Step1.
  # lmersamp <- lmer(yij~x + (1|area), data = datsamp) #f(yij|i=1,j=1)
  # sig2ehatsamp <- summary(lmersamp)$sigma^2
  
  if(isSingular2(lmersamp)){
    #Adjusted Likelihood
    sig2ehatsamp <- summary(lmersamp)$sigma^2
    sig2uhatsamp <- optim(par=VarCorr(lmersamp)[[1]][1],fn=P.LikHood,datsamp=datsamp,method="Brent",lower=0,upper=200,lmersamp=lmersamp)$par
  }else{
    ##### Implement EBUP simulation method:
    sig2ehatsamp <- summary(lmersamp)$sigma^2
    sig2uhatsamp <- VarCorr(lmersamp)[[1]][1]
  }
  
  # method 1 linear regression with estimated ui
  ui <- ranef(lmersamp)[[1]][,1];  names(ui) <- rownames(ranef(lmersamp)[[1]])
  datsamp$est.ui <- ui[as.factor(datsamp$area)];init.coef <- lm(log(1/probi)~est.ui,data=datsamp)
  Init <- c(init.coef$coef,summary(init.coef)$sigma)
  
  # Gaussian Hermite approximation.
  GH_rule <- gaussHermiteData(M)
  U = GH_rule$x
  W = GH_rule$w
  log_pi = log(tapply(datsamp$probi, datsamp$area,function(.){sample(.,size=1)}) )
  
  toptim <- optim(Init, function(.){Fct(., W, U, log_pi, uhatis, vcondus)})
  
  # toptim$convergence
  est_par = toptim$par
  est_par[3] = est_par[3]^2
  names(est_par) = paste0("par",1:3)
  
  return(list(par = est_par,
              convergence = toptim$convergence))
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
  
  
  #################################################################################
  
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
    
    # #############################Coefarea parameters  08022021 editted################################################################
    # coefarea.initial <- Coef_Area_GH(datsamp,lmersamp,M=100,uhatis, vcondus)$par#
    # ###Closed formula
    # coefarea <- optim(c(coefarea.initial[1:2],sqrt(coefarea.initial[3])),parweight, uhatis= uhatis, vcondus = vcondus, datsamp = datsamp)$par
    # coefarea[3] <- (coefarea[3])^2
    # # 
    #######################################################################################################################
    
    
    theta.hat =  list(betahat = betahatsamp
                      , sig2uhat = sig2uhatsamp
                      , sig2ehat = sig2ehatsamp
                      , xbar = xbaris
                      , ybar = ybaris
                      , gammahat = gammahatis 
                      , coefunit_12 = coefunit[1:2]
                      # ,  coefarea =  coefarea
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



gen.vcov.jack = function(datsamp, par_method){
  
  sampnum = length(unique(datsamp$area))
  samparea = unique(datsamp$area)   # updated 08192021!!!!!!!
  
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





parweight = function(w.coef, uhatis,vcondus,datsamp){         # Use the EBUP document page 5. f_{s}(u_i|D_s,I_i = 1)                     ######## added for the nonsampled areas
  iprob =  tapply(datsamp$probi, datsamp$area,function(.){.[1]})
  weight = 1/iprob;  names(weight)=names(uhatis);
  gam1 = w.coef[1];gam2 = w.coef[2];tau = w.coef[3]
  
  new.var = (tau/gam2)^2 +vcondus
  new.mean = ((log(weight)-gam1)/gam2 -uhatis)^2
  
  intral = -log(abs(gam2))-log(2*pi*new.var)/2-new.mean/(2*new.var)
  # intral = igam2*(1/sqrt(2*pi*new.var))*exp(-(mu1/(2*new.var)))
  
  likfc = -sum(intral)
} 


gen.vcov.jack.area = function(datsamp,uhatis,vcondus){
  samparea= unique(datsamp$area)
  sampnum = length(unique(datsamp$area))
  
  par.mat = lapply(1:sampnum, function(b){
    datsamp_complement = datsamp[datsamp$area!=samparea[b],]
    lmersamp_complement <- lmer(yij~x + (1|area), data = datsamp_complement,REML=FALSE ) 
    
    ###Coefarea parameters  08022021 editted
    est.ui_comp <- ranef(lmersamp_complement)$area; names(est.ui_comp) = "est.ui"
    new.dat_comp = cbind(probi= tapply(datsamp_complement$probi, datsamp_complement$area,function(.){.[1]}) , est.ui_comp)
    lm.ini_comp = lm(log(1/probi)~est.ui, data = new.dat_comp)
    coefarea.initial_comp = c(lm.ini_comp$coefficients, var(log(1/new.dat_comp$probi)))
    # coefarea.initial <- Coef_Area_GH(datsamp,lmersamp,M=100,uhatis,vcondus)$par#
    
    ###Closed formula
    coefarea_comp <- optim(c(coefarea.initial_comp[1:2],sqrt(coefarea.initial_comp[3])),parweight, uhatis= uhatis[-b], vcondus = vcondus[-b], datsamp = datsamp_complement)$par
    coefarea_comp[3] <- (coefarea_comp[3])^2
    coefarea_comp  
  }) %>% do.call("rbind",.)
  par.mean = apply(par.mat,2,mean)
  par.sum = lapply(1:sampnum, function(i){
    
    diff = as.matrix( par.mat[i,]-par.mean)
    diff.mat = diff%*%t(diff)
  })
  return((sampnum-1)* Reduce("+",par.sum)/sampnum)
}


dat.gen = function(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, par_method = "Ori",stratindarea, stratindpop){
  
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
  
  #### estimated unit weight model parameters  ### CHANGED PART IT WORKS WELL for JACKKNIFE 08192021.R
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
    coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 500000))$par #11292021 increase maxit number
  }
  # coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 10000))$par
  
  #update 02092021
  vcov.jack = gen.vcov.jack(datsamp, par_method)
  
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
  
  
  #################################################################
  ###Coefarea parameters  08022021 editted
  est.ui <- ranef(lmersamp)$area; names(est.ui) = "est.ui"
  new.dat = cbind(probi= tapply(datsamp$probi, datsamp$area,function(.){.[1]}) , est.ui)
  lm.ini = lm(log(1/probi)~est.ui, data = new.dat)
  coefarea.initial = c(lm.ini$coefficients, var(log(1/new.dat$probi)))
  # coefarea.initial <- Coef_Area_GH(datsamp,lmersamp,M=100,uhatis,vcondus)$par#
  
  ###Closed formula
  coefarea <- optim(c(coefarea.initial[1:2],sqrt(coefarea.initial[3])),parweight, uhatis= uhatis, vcondus = vcondus, datsamp = datsamp)$par
  coefarea[3] <- (coefarea[3])^2
  # 
  
  vcov.jack.area = gen.vcov.jack.area(datsamp,uhatis,vcondus)
  
  # 
  # ui <- ranef(lmersamp)[[1]][,1];  names(ui) <- rownames(ranef(lmersamp)[[1]])
  # datsamp$est.ui <- ui[as.character(datsamp$area)];init.coef <- lm(log(1/probi)~est.ui,data=datsamp)
  # # method 1 linear regression with estimated ui
  # Init <- c(init.coef$coef,summary(init.coef)$sigma)
  # 
  # coefarea <- optim(c(Init[1:2],sqrt(Init[3])),parweight, uhatis= uhatis, vcondus = vcondus, datsamp = datsamp)$par
  # coefarea[3] <- (coefarea[3])^2
  # coefarea
  if(par_method == "Ori"){
    par.mu = c(betahatsamp, sig2uhatsamp, sig2ehatsamp, coefunit[1], coefunit[2], coefarea)  ### coefarea added 08022021
  }else if (par_method == "PS"){
    par.mu = c(betahatsamp, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0, coefarea)
  }else if (par_method == "AW"){
    par.mu = c(betahatsamp.aw, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0, coefarea)
  }else if (par_method == "W"){
    par.mu = c(betahatsamp.w, sig2uhatsamp, sig2ehatsamp, coefunit[1], 0, coefarea)
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
       , coefarea = coefarea
       , vcov.jack.area = vcov.jack.area
  )

}


### Added 07302022 for Dr. Molina's MSE estimator:


par.mu.est = function(datsamp){
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
  
  #### estimated unit weight model parameters  ### CHANGED PART IT WORKS WELL for JACKKNIFE 08192021.R
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
    coefunit <- optim(lmpireg$coef, sseexppi, ypi = 1/datsamp$probij ,y = datsamp$yij, x = datsamp$x, areaspi= datsamp$area,control=list(maxit = 500000))$par #11292021 increase maxit number
  }
  ###Coefarea parameters  08022021 editted
  est.ui <- ranef(lmersamp)$area; names(est.ui) = "est.ui"
  new.dat = cbind(probi= tapply(datsamp$probi, datsamp$area,function(.){.[1]}) , est.ui)
  lm.ini = lm(log(1/probi)~est.ui, data = new.dat)
  coefarea.initial = c(lm.ini$coefficients, var(log(1/new.dat$probi)))
  # coefarea.initial <- Coef_Area_GH(datsamp,lmersamp,M=100,uhatis,vcondus)$par#
  
  ###Closed formula
  coefarea <- optim(c(coefarea.initial[1:2],sqrt(coefarea.initial[3])),parweight, uhatis= uhatis, vcondus = vcondus, datsamp = datsamp)$par
  coefarea[3] <- (coefarea[3])^2
  # 
  
  par.mu = c(betahatsamp, sig2uhatsamp, sig2ehatsamp, coefunit[1], coefunit[2], coefarea) 
  
  
  par.mu
  
}
######################################
# Fct <- function(coef,W,U,log_pi,sig2uhatsamp, uhatis, vcondus){
#   
#   par1 = coef[1]
#   par2 = coef[2]
#   tau = abs(coef[3]); itau = 1/tau
#   sigma = sqrt(sig2uhatsamp)
#   
#   Mat = W*itau*exp( -((-log_pi-par1-par2*sqrt(2)*sigma*U)/(sqrt(2)*tau))^2)
#   
#   Neg_ll = - sum(log(apply(Mat, 1, sum)))
#   
#   return(Neg_ll)
# }
# # Gaussian Hermite method is used.
# Coef_Area_GH = function(datsamp,lmersamp, M = 100){
#   # Step1.
#   # lmersamp <- lmer(yij~x + (1|area), data = datsamp) #f(yij|i=1,j=1)
#   # sig2ehatsamp <- summary(lmersamp)$sigma^2
#   
#   if(isSingular2(lmersamp)){
#     #Adjusted Likelihood
#     sig2ehatsamp <- summary(lmersamp)$sigma^2
#     sig2uhatsamp <- optim(par=VarCorr(lmersamp)[[1]][1],fn=P.LikHood,datsamp=datsamp,method="Brent",lower=0,upper=200,lmersamp=lmersamp)$par
#   }else{
#     ##### Implement EBUP simulation method:
#     sig2ehatsamp <- summary(lmersamp)$sigma^2
#     sig2uhatsamp <- VarCorr(lmersamp)[[1]][1]
#   }
#   
#   # method 1 linear regression with estimated ui
#   ui <- ranef(lmersamp)[[1]][,1];  names(ui) <- rownames(ranef(lmersamp)[[1]])
#   datsamp$est.ui <- ui[as.factor(datsamp$area)];init.coef <- lm(log(1/probi)~est.ui,data=datsamp)
#   Init <- c(init.coef$coef,summary(init.coef)$sigma)
#   
#   # Gaussian Hermite approximation.
#   GH_rule <- gaussHermiteData(M)
#   length_data = nrow(datsamp)
#   U = kronecker(t(GH_rule$x), rep(1,length_data))
#   W = kronecker(t(GH_rule$w), rep(1,length_data))
#   log_pi = kronecker(log(datsamp$probi), t(rep(1,M)))
#   
#   toptim <- optim(Init, function(.){Fct(., W, U, log_pi, sig2uhatsamp)})
#   
#   # toptim$convergence
#   est_par = toptim$par
#   est_par[3] = est_par[3]^2
#   names(est_par) = paste0("par",1:3)
#   
#   return(list(par = est_par,
#               convergence = toptim$convergence))
# }
