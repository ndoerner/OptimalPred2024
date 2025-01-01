
# This code is modified from fns_MSE_estimation_0729.R
# pop.wt.bt =  rep(kappai, Nis)*exp(par.mu[5]*datsamparea$x+par.mu[6]*pop.yij.bt) #07182021 changed  in EBP fully
# EBP_fully two version
# MSE_PS_1iter editted
# MSE fully will be not used in final scenarios!!!!!!!!!! noting editted


# This code is modified from fns_MSE_estimation_nonsamp_0805.R
## cHANGED THETA.B TO GENERATE POSITIVE sig2uhat.b : Delta Method.

# Updated from 10262021/fns_MSE_estimation_nonsamp_0819.R
##

source("fns_predictors_nonsamp_0819.R")
source("fns_datasets_nonsamp_0819.R")

########################################################################################################
# 2021.06.22 update : make functions for making bc MSE estimator

EBP_nonsamp_infor <- function(BL, par.mu,vcov.jack,Iareas,datnonsamparea, vcov.jack.area,xbaris, ybaris, nis){
  
  
    ## editted 08192021 Delta Method for Jackknife Variance Estimation##
    par.mu.prim = c(par.mu[1],par.mu[2],log(par.mu[3]), log(par.mu[4]),par.mu[5],par.mu[6])
    par.mu.jacobian = diag(c(1, 1, 1/par.mu[3], 1/par.mu[4], 1, 1)) 
    vcov.jack.prim = t(par.mu.jacobian) %*% vcov.jack %*% (par.mu.jacobian)
    theta.b <- mvrnorm(n = 1, mu = par.mu.prim , Sigma = vcov.jack.prim)
    theta.b[3] = exp(theta.b[3])
    theta.b[4] = exp(theta.b[4])
    
    # betahatsamp <- fixef(lmersamp)
    # nis <- tapply(datsamp$area, datsamp$area, length ) # should be length!!!!!!!!!!! not sum!!!!!!!!!!!
    # n <-sum(nis)
    # ybaris <- tapply(datsamp$yij, datsamp$area, mean)
    # xbaris <- tapply(datsamp$x, datsamp$area, mean)
    gammahatis.b <- theta.b[3]/(theta.b[3] + theta.b[4]/nis)
    uhatis.b <- gammahatis.b*(ybaris - theta.b[1] - theta.b[3][2]*xbaris)
    vcondus.b <- gammahatis.b*theta.b[4]/nis
    
    # par.mu.area = c(par.mu[7], par.mu[8], log(par.mu[9]))
    # par.mu.area.jacobian = diag(c(1,1,1/par.mu[9]))
    # vcov.jack.area.prim = t(par.mu.area.jacobian) %*% vcov.jack.area %*% (par.mu.area.jacobian)
    coefarea.b <- mvrnorm(n=1, mu = par.mu[7:8], Sigma = vcov.jack.area[1:2,1:2])
    coefarea.b <- c(coefarea.b, par.mu[9])
    # coefarea.b[3] = exp(coefarea.b[3])
    # coefarea.b
    # 
    yns <- replicate(BL, gendata_nonsamp_Inversion(Iareas,coefarea.b, theta.b[3],theta.b[1:2],datnonsamparea, theta.b[4],theta.b[6])) #nrow(datsamparea)*BL
    # yns[datsamparea$Iunits==1] <- datsamparea$yij[datsamparea$Iunits==1]
    
    gc()  
  
  return(t(yns))
}


Informative_nonsamp_BC <- function(BL, par.mu,vcov.jack,Iareas,datnonsamparea, vcov.jack.area,xbaris, ybaris, nis, B){

  # eg . Two same codes :  
  # lapply(1:3, function(.){matrix(1:4, nrow = 2)}) %>% sapply(., function(res){res}); sapply(1:3, function(.){matrix(1:4, nrow = 2)});
  
  yNs <- sapply(1:B, function(i){
    # T.vec = NaN 
    # while(!all(!is.na(T.vec))){  #muted 08192021 thanks to delta method!
    T.vec = EBP_nonsamp_infor(BL, par.mu,vcov.jack,Iareas,datnonsamparea, vcov.jack.area,xbaris, ybaris, nis) 
    # }
    return(T.vec)
  })
  rm(datnonsamparea)
  gc()
  return(yNs)
  
}

MSE_nonsamp_BC <- function(BL, par.mu,vcov.jack,Iareas,datnonsamparea,B, hlist,EBP, vcov.jack.area,xbaris, ybaris, nis){
  nh = length(hlist)
  T.par = paste0("h",1:nh) #changed 08192021
  
  yetaboot <- Informative_nonsamp_BC(BL, par.mu,vcov.jack,Iareas,datnonsamparea, vcov.jack.area,xbaris, ybaris, nis,B)
  
  colnames(yetaboot) <-1:B ##06292021
  
  gc()
  yetaboot <- yetaboot %>% cbind.data.frame(
    rep = rep(1:BL, times = nrow(datnonsamparea)),
    area = rep(datnonsamparea$area, each = BL)
  )
  gc()
  
  theta_rep_boot <- lapply(1:B, function(i){
    
      yetaboot %>% dplyr::select(i, area, rep) %>% 
      # mutate_at(as.character(b), hlist) %>%
      setNames(c("theta", "area", "rep")) %>%
      group_by(area, rep) %>%  ## order has been changed
      mutate(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), h5= hlist[["h5"]](theta), h6= hlist[["h6"]](theta), h7= hlist[["h7"]](theta)) %>%  # h8 added 08192021
      # mutate(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), h5= hlist[["h5"]](theta), h6= hlist[["h6"]](theta), h7= hlist[["h7"]](exp(theta))) %>%  #gini coef ; use exp(yij) editted 07302021
      group_by(area, rep) %>%
      summarise_at(vars(contains("h")), mean) %>% 
      mutate(boot = i) %>% 
      return(.)
    
  }) %>% do.call("rbind",.)
  
  rm(yetaboot) ## 07012021 editted
  gc()
  
  
  M1hat_boot <- theta_rep_boot %>% dplyr::select(-rep) %>% 
    group_by( boot, area)  %>% summarise_all(var)  %>% group_by(area) %>% summarise_at(vars(contains("h")), mean) %>% 
    dplyr::select(T.par)%>% as.matrix(.)  # return 120*6 matrix
  
  
  ## Boot M2 ####
  estpssampa_boot <- theta_rep_boot %>% dplyr::select(-rep) %>% group_by(area,boot) %>% 
    summarise_all(mean) %>% arrange(boot,area) %>% ungroup()
  gc()
  
  # rm(theta_rep_boot) ## 07012021 editted
  # gc()
  
  M2hat_boot <- sapply(1:nh, function(hi){
    colname <- paste0("h", hi)
    tapply((unlist(estpssampa_boot[, colname]) - (EBP[, colname]))^2, estpssampa_boot$area, mean)
  })
  
  colnames(M2hat_boot) = paste0("h",1:nh)
  
  return(list( M1hat_bt = M1hat_boot,
               M2hat_bt = M2hat_boot,
               theta_rep_boot = theta_rep_boot # added 06092022 for CondCI
               ))
  
}


EBP_ML <- function(Iareas,datsamp, datnonsamparea, par_method = "Ori", BL,par.mu, hlist){
  
  nh = length(hlist)
  T.par = paste0("h",1:nh) #CHanged
  Nis = tapply(datnonsamparea$yij, datnonsamparea$area, length)
  
  n.nonsamparea = length(unique(datnonsamparea$area))
  betahat = c(par.mu[1], par.mu[2]);sig2uhat = par.mu[3];sig2ehat = par.mu[4];ahat = par.mu[5];bhat = par.mu[6]
  lambda1hat = par.mu[7];lambda2hat = par.mu[8];tau2hat = par.mu[9]
  
  ## Generate Bootstrap Nonsampled Areas:
  # pop.yij.b <- gendata_nonsamp_Inversion(Iareas,coefarea,sig2uhatsamp, betahatsamp,datnonsamparea, sig2ehatsamp,b)
  pop.yij.b <- gendata_nonsamp_Inversion(Iareas,par.mu[7:9], par.mu[3], par.mu[1:2], datnonsamparea, par.mu[4], par.mu[6])
  datnonsamparea.b = datnonsamparea
  datnonsamparea.b$yij = pop.yij.b
  
  rm(datnonsamparea);gc()
  
  ## Step2 : Bootstrap version True Par:
  
  theta.b = sapply(1:nh, function(i){tapply(datnonsamparea.b$yij, datnonsamparea.b$area, hlist[[i]])}) %>%
    as.data.frame(.) %>% setNames(paste0("h",1:nh)) 
  
  ## Step 3 : Generate Bootstrap sample (pop.yij.samp.b) \& Model parameter Estimates
  n.samparea = length(unique(datsamp$area))
  
  uhat.b <- rnorm(n.samparea, mean = 0, sd = sqrt(sig2uhat)) %>% setNames( names(Iareas[Iareas==1]))
  pop.yij.samp.b <- bhat*sig2ehat + betahat[1] + betahat[2]*datsamp$x + uhat.b[as.character(datsamp$area)] + rnorm(nrow(datsamp), mean = 0, sd = sqrt(sig2ehat))
  
  datsamp.b = datsamp
  datsamp.b$yij = pop.yij.samp.b
  
  par.mu.b = par.mu.est(datsamp.b)
  
  theta.hat.b = EBP_nonsamp_MC(BL, Iareas, par.mu.b[7:9], par.mu.b[3], par.mu.b[1:2], datnonsamparea.b, par.mu.b[4], par.mu.b[5:6], par_method = "Ori", hlist)$estpssampa%>% as.matrix(.) %>% .[,T.par]
  
  
  MSE_bt = ((theta.hat.b - theta.b)^2) %>% as.matrix(.)
  apply(MSE_bt,2,mean)
  rm(theta.hat.b, theta.b);gc()
  return(MSE_bt)
}




MSE_ML <- function(Iareas,datsamp, datnonsamparea, par_method = "Ori", BL,par.mu, hlist,B){
  
  MSE_ys <- lapply(1:B, function(i){
    
    EBP_ML(Iareas,datsamp, datnonsamparea, par_method = "Ori", BL,par.mu, hlist)
    
  })%>% Reduce("+",.)
  
  MSE_ml = MSE_ys/B
  rm(datsamp,datnonsamparea)
  gc()
  return(MSE_ml)
  
}





# #######################################################################################################################
# 
# EBP_jack = function(i,datsamp, datsamparea, samparea, par_method="Ori",BL, nis, xbaris, ybaris){
#   
#   pars <- gen.est.par(datsamp[datsamp$area!=samparea[i],], par_method)
#   pars.b <- c(pars$betahat,pars$sig2uhat, pars$sig2ehat, pars$coefunit_12)
#   betahatsamp.b <- pars.b[1:2]; sig2uhatsamp.b <- pars.b[3]; sig2ehatsamp.b<- pars.b[4]; bs<-pars.b[6]
#   gammahatis.b <- sig2uhatsamp.b/(sig2uhatsamp.b + sig2ehatsamp.b/nis)
#   uhatis.b <- gammahatis.b*(ybaris - betahatsamp.b[1] - betahatsamp.b[2]*xbaris)
#   vcondus.b <- gammahatis.b*sig2ehatsamp.b/nis
#   yNs <- replicate(BL, gendata(uhatis.b,vcondus.b,betahatsamp.b,datsamparea,sig2ehatsamp.b,bs))
# 
#   yNs[datsamparea$Iunits==1]<-datsamparea$yij[datsamparea$Iunits==1]
#   return(t(yNs))
# }
# 
# Informative_JACK <- function(datsamp, datsamparea,samparea, par_method = "Ori",BL, nis, xbaris, ybaris, sampnum){
#   
#   
#   yNs = sapply(1:sampnum, function(i){
#     EBP_jack(i,datsamp, datsamparea,samparea, par_method="Ori",BL, nis, xbaris, ybaris)
#   })
#   rm(datsamp, datsamparea)
#   gc()
#   return(yNs)
# 
# }
# 
# MSE_Jack <- function(datsamp, datsamparea, samparea, par_method="Ori",BL, nis, xbaris, ybaris, hlist, EBP, sampnum){
#   nh = length(hlist)
#   yetajack <- Informative_JACK(datsamp, datsamparea,samparea, par_method = "Ori",BL, nis, xbaris, ybaris, sampnum)
#   colnames(yetajack) <-1:sampnum
#   
#   yetajack <- yetajack %>% cbind.data.frame(
#     area = rep(datsamparea$area, each = BL),
#     rep = rep(1:BL, times = nrow(datsamparea))
#   )
#   gc()
#   
#   theta_rep_jack <- lapply(1:sampnum, function(i){
#     
#     yetajack %>% dplyr::select(i, area, rep) %>% 
#       # mutate_at(as.character(b), hlist) %>%
#       group_by(area, rep) %>%  ## order has been changed
#       setNames(c("theta", "area", "rep")) %>%    
#       mutate(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), h5= hlist[["h5"]](theta), h6= hlist[["h6"]](theta), h7= hlist[["h7"]](exp(theta))) %>%  #gini coef ; use exp(yij) editted 07302021
#       group_by(area, rep) %>%
#       summarise_at(vars(contains("h")), mean) %>% 
#       mutate(boot = samparea[i]) %>%
#       return(.)
#     
#   }) %>% do.call("rbind", .)
#   
#   rm(yetajack) ## 07012021 editted
#   gc()
#   
#   M1hat_jack <- theta_rep_jack %>% dplyr::select(-rep) %>% 
#     group_by(area,boot) %>% summarise_all(var) 
#   
#   M1hat.jk.Rao = sapply(1:nh, function(i){
#     colnames = paste0("h",i)
#     # bias = (sampnum-1)*tapply(pull(M1hat_boot[,colnames])-rep((M1.bt[,colnames]),each =sampnum), M1hat_boot$area,mean)
#     bias = (sampnum-1)* tapply(pull(M1hat_jack[,colnames]), M1hat_jack$area,mean)
#     bias
#   }) %>% as.data.frame(.) %>% setNames(paste0("h",1:nh))
#   gc()
#   
#   M1hat.jk.Lohr = sapply(1:nh, function(i){
#     colnames=paste0("h",i)
#     M1hat_i = M1hat_jack%>%dplyr::select(area,boot,colnames)%>%filter(.,boot!=area)
#     bias = tapply(pull(M1hat_i[,colnames]), M1hat_i$area,sum)
#     # bias = tapply(pull(M1hat_i[,colnames])-rep(pull(M1hat[,colnames]),each =(sampnum-1)), M1hat_i$area,sum)
#   })%>%as.data.frame(.)%>%setNames(paste0("h",1:nh))
#   gc()
#   
#   ## bootstrap EBUP ####
#   estpssampa_jack <- theta_rep_jack %>% dplyr::select(-rep) %>% group_by(area, boot) %>%
#     summarise_all(mean) %>% arrange(boot, area) %>% ungroup()
#   rm(theta_rep_jack) ## 07012021 editted
#   gc()
#   
#   
#   ## Boot M2 ####
#   M2hat_jack <- sapply(1:nh, function(hi){
#     colname <- paste0("h", hi)
#     (sampnum-1)*tapply((unlist(estpssampa_jack[, colname]) - (EBP[, colname]))^2, estpssampa_jack$area, mean)
#   })
#   colnames(M2hat_jack)=paste0("h",1:nh)
#   
#   M2hat.jack <- as.data.frame(M2hat_jack) %>%mutate(area=samparea) %>%.[,c(paste0("h",1:nh))]
#   
#   gc()
#   return(list(M1hat.jk.Rao = M1hat.jk.Rao,
#               M1hat.jk.Lohr = M1hat.jk.Lohr,
#               M2hat.jk = M2hat.jack))
# }
# 
# ####################################################################################################################################
# #################################Version1
# 
# EBP_fully <- function(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist){
#   nh = length(hlist)
#   T.par = c("h1","h2","h3","h4","h5","h6","h7")
#   Nis = tapply(datsamparea$yij, datsamparea$area, length)
#   ## Genereate Bootsrap Parameter estimation
#   cond.moments = cond.mu.v(par.mu[1:2], par.mu[3], par.mu[4], xbaris, ybaris, gammahatis)
#   T.uhatis = cond.moments[[1]]; T.vcondus = cond.moments[[2]]
#   pop.yij.bt = gendata(T.uhatis , T.vcondus,  par.mu[1:2], datsamparea , par.mu[4], b=par.mu[6])
#   
#   kappai = coefunit[2+1:length(samparea)]
#   hat.a = par.mu[5]
#   pop.wt.bt =  rep(kappai, Nis)*exp(par.mu[5]*datsamparea$x+par.mu[6]*datsamparea$yij)
#   
#   
#   ## Estimate 
#   datsamparea2 = datsamparea
#   datsamparea2$yij = pop.yij.bt
#   pop.theta = sapply(hlist, function(h) tapply((datsamparea2$yij), datsamparea2$area, h))%>%  #0713 changed for the new parameters
#     as.data.frame(.) %>% setNames(paste0("h",1:nh)) 
#   
#   datsamparea2$probij = 1/pop.wt.bt
#   datsamp2 = datsamparea2[datsamparea2$Iunits ==1,]
#   par.b =  gen.est.par(datsamp2, par_method = "Ori")
#   coefunit2 = par.b$coefunit
#   
#   hat.theta = EBP_MC(BL, par.b$betahat, par.b$sig2uhat, par.b$sig2ehat, par.b$xbar, par.b$ybar, par.b$gammahat, datsamparea2, coefunit2, par_method = "Ori", hlist)$estpssampa %>%
#     as.matrix(.) %>% .[,T.par]
#   
#   
#   MSE_bt = ((hat.theta - pop.theta)^2) %>% as.matrix(.)
#   rm(hat.theta) ## 07012021 editted
#   gc()
#   return(MSE_bt)
# }
# 
# 
# 
# 
# MSE_Fully <- function(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist,B){
#   
#   MSE_ys <- sapply(1:B, function(i){
#     
#     EBP_fully(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist)
#     
#   })
#   rm(datsamp,datsamparea)
#   gc()
#   return(MSE_ys)
#   
# }
# 
# 
# 
# 
# 
# ####################################################################################################################################
# ################################# Version2
# 
# 
# EBP_fully2 <- function(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist){
#   nh = length(hlist)
#   T.par = c("h1","h2","h3","h4","h5","h6","h7")
#   Nis = tapply(datsamparea$yij, datsamparea$area, length)
#   ## Genereate Bootsrap Parameter estimation
#   cond.moments = cond.mu.v(par.mu[1:2], par.mu[3], par.mu[4], xbaris, ybaris, gammahatis)
#   T.uhatis = cond.moments[[1]]; T.vcondus = cond.moments[[2]]
#   pop.yij.bt = gendata(T.uhatis , T.vcondus,  par.mu[1:2], datsamparea , par.mu[4], b=par.mu[6])
#   
#   kappai = coefunit[2+1:length(samparea)]
#   hat.a = par.mu[5]
#   pop.wt.bt =  rep(kappai, Nis)*exp(par.mu[5]*datsamparea$x+par.mu[6]*pop.yij.bt) ##############################07182021 changed from fns_MSE_estimation_0715.R
#   
#   
#   ## Estimate 
#   datsamparea2 = datsamparea
#   datsamparea2$yij = pop.yij.bt
#   pop.theta = sapply(hlist, function(h) tapply((datsamparea2$yij), datsamparea2$area, h))%>%  #0713 changed for the new parameters
#     as.data.frame(.) %>% setNames(paste0("h",1:nh)) 
#   
#   datsamparea2$probij = 1/pop.wt.bt
#   datsamp2 = datsamparea2[datsamparea2$Iunits ==1,]
#   par.b =  gen.est.par(datsamp2, par_method = "Ori")
#   coefunit2 = par.b$coefunit_12  # MISTAKE
#   
#   hat.theta = EBP_MC(BL, par.b$betahat, par.b$sig2uhat, par.b$sig2ehat, par.b$xbar, par.b$ybar, par.b$gammahat, datsamparea2, coefunit2, par_method = "Ori", hlist)$estpssampa %>%
#     as.matrix(.) %>% .[,T.par]
#   
# 
#   MSE_bt = ((hat.theta - pop.theta)^2) %>% as.matrix(.)
#   rm(hat.theta) ## 07012021 editted
#   gc()
#   return(MSE_bt)
# }
# 
# 
# 
# 
# MSE_Fully2 <- function(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist,B){
#   
#   MSE_ys <- sapply(1:B, function(i){
#     
#     EBP_fully2(datsamp, datsamparea, par_method = "Ori", BL, par.mu, xbaris, ybaris, gammahatis, coefunit, samparea,hlist)
#     
#   })
#   rm(datsamp,datsamparea)
#   gc()
#   return(MSE_ys)
#   
# }
# 


#################################################################################################################
# MSE_PS_1iter changed 08032021
MSE_PS_1iter = function(datsamp, samparea, par.mu, datnonsamparea,par_method = "Ori" , Iareas, coefunit){
  
  ## Step1- step2
  datsamp2 = datsamp
  u_samp_b = rnorm(length(samparea), mean = 0, sd = sqrt(par.mu[3]))
  names(u_samp_b) = as.character(samparea)
  nis = tapply(datsamp$yij, datsamp$area, length)
  samp.yij.bt = par.mu[1] + par.mu[2]*datsamp$x + rep(u_samp_b, nis) + rnorm(sum(nis), mean = 0, sd = sqrt(par.mu[4])) ###########sum(nis)
  
  datsamp2$yij = samp.yij.bt ##################changed
  
  kappai = coefunit[2+1:length(samparea)]
  hat.a = par.mu[5]
  samp.wt.bt =  rep(kappai, nis)*exp(par.mu[5]*datsamp2$x + par.mu[6]*datsamp2$yij)
  
  datsamp2$probij = 1/samp.wt.bt  ##################changed
  
  ## Step3
  ## new parameter estimates
  par.b =  gen.est.par(datsamp2, par_method = "Ori")
  uhatis.b = par.b$gammahat*(par.b$ybar - par.b$betahat[1]-par.b$betahat[2]*par.b$xbar) 
  nonsamparea.b = EBP_nonsamp_PS(datsamp, datnonsamparea, Iareas, par.b$betahat,  par.b$coefunit_12[2], par.b$sig2uhat, uhatis.b)
  
  
  rm(datsamp, datnonsamparea)
  gc()
  return(nonsamparea.b)
}


MSE_PS = function(datsamp, samparea, par.mu, datnonsamparea,par_method = "Ori" , Iareas, coefunit, B){
  Boots_PS_B = sapply(1:B, function(i){
    MSE_PS_1iter(datsamp, samparea, par.mu, datnonsamparea,par_method = "Ori" , Iareas, coefunit)
  })
  Nis = tapply(datnonsamparea$yij, datnonsamparea$area, length)
  MSE.PS1 = apply(Boots_PS_B,1,var)*(B-1)/B
  MSE.PS2 = par.mu[4]/Nis
  
  gami = tapply(datsamp$yij-par.mu[1]-par.mu[2]*datsamp$x, datsamp$area, mean)
  W = tapply(datsamp$probi,datsamp$area,function(.){sample(.,size = 1)}); wi=1/W; names(wi)<-names(Iareas[Iareas==1])
  Swi_1 = sum(wi-1)
  nis = tapply(datsamp$yij, datsamp$area, length)
  MSE.PS3 = sum(((wi-1)*(gami-sum(((wi-1)*gami/Swi_1)))^2)/Swi_1)-sum(par.mu[4]/nis)/length(nis)
  
  MSE_PS_h1 = MSE.PS1 + MSE.PS2 + MSE.PS3
  
  nonsamparea = Iareas[Iareas ==0]
  names(MSE_PS_h1) = names(nonsamparea)
  rm(Boots_PS_B) ## 07022021 editted
  gc()
  
  return(MSE_PS_h1)
}



############all sampled area

com.vec= function(Iareas,vec){
  Noareas = Iareas[Iareas==0]; Novec = rep(NA,length(Noareas)); names(Novec)=names(Noareas)
  names(vec) = names(Iareas[Iareas==1])
  allareamean = c(vec,Novec)
  allareamean =allareamean[match(names(Iareas),names(allareamean))]
  allareamean
}



