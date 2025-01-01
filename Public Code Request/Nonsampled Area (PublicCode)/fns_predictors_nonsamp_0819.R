# This code is modified from functions_MSE.R
# This code is modified from fns_predicots_0623.R : parameters are changed!!!! ex. Gini coefficient
# This code is modified from fns_predicots_0715.R : use the  "exp(yij)" for Gini coefficient
# This code is modified from fns_predicots_0805.R : Inversion of F_{c}(u_i)
# This code is modified from fns_predicots_0805.R : changed in EBP_MC for including function "h8" : gini(abs(y)) 08192021
# This code is modified from fns_predicots_0819.R : changed in FsUi becasue use wrong assumption on f_{s}(u_i)

## EBP calculation by Methods:
cond.mu.v = function(betahat, sig2uhatsamp, sig2ehatsamp, xbaris.w, ybaris.w, gammahatis.w ){
  uhatis.w <- gammahatis.w*(ybaris.w - betahat[1] - betahat[2]*xbaris.w)
  vcondus.w <- (1-gammahatis.w)*sig2uhatsamp
  return(list(uhatis.w = uhatis.w
              , vcondus.w = vcondus.w))
  
}

######################################################################################################################
############       Chaged from this part :: do not use "gendata"::generate f_{c}(u_i)         #########################
######################################################################################################################

simynonsamp <-function(nonareaui,betahatsamp,datnonsamparea,sig2ehatsamp,b){
 
  ycondgen <- b*sig2ehatsamp + betahatsamp[1] + betahatsamp[2]*datnonsamparea$x + nonareaui[as.character(datnonsamparea$area)]+ rnorm(nrow(datnonsamparea), mean = 0, sd = sqrt(sig2ehatsamp))
  return(ycondgen)
  
}


FsUi = function(ui,coefarea,sig2uhatsamp){              ####changed 09162021 FcUi is the correct representation:(
  gam1 = coefarea[1];gam2 = coefarea[2];tau2=coefarea[3]
  Denom= exp(gam1+tau2/2+(sig2uhatsamp*gam2^2)/2)-1
  mult = exp(gam1+tau2/2+(sig2uhatsamp*gam2^2)/2)
  Numer = mult*pnorm(ui,mean=(sig2uhatsamp*gam2),sd= sqrt(sig2uhatsamp))-pnorm(ui,mean=0,sd=sqrt(sig2uhatsamp))## changed mean!! from -(sig2uhatsamp*gam2) to (sig2uhatsamp*gam2)
  return(Numer/Denom)
}


gendata_nonsamp_Inversion<-function(Iareas,coefarea,sig2uhatsamp, betahatsamp,datnonsamparea, sig2ehatsamp,b){
  
  nonarea = Iareas[Iareas!=1]
  Fs <- function(x) { FsUi(x,coefarea=coefarea,sig2uhatsamp=sig2uhatsamp)}
  Fs <- Vectorize(Fs)
  F.inv <- function(y){uniroot(function(x){Fs(x)-y},interval=c(-1000000,1000000))$root}
  F.inv <- Vectorize(F.inv)
  X <- runif(length(nonarea),0,1)   # random sample from U[0,1]
  nonareaui <- F.inv(X)
  names(nonareaui) <- names(nonarea)
  
  ###document model
  ypsimnonsamparea<- simynonsamp(nonareaui,betahatsamp,datnonsamparea,sig2ehatsamp,b)
  
  return( ypsimnonsamparea = ypsimnonsamparea)
}


# gendata <- function(uhatis,vcondus,betahatsamp,datsamparea,sig2ehatsamp,b){
#   ucondgen <- rnorm(length(uhatis), mean = uhatis, sd = sqrt(vcondus))
#   names(ucondgen) <- names(uhatis)
#   ycondgen <- b*sig2ehatsamp + betahatsamp[1] + betahatsamp[2]*datsamparea$x + ucondgen[as.character(datsamparea$area)] + rnorm(nrow(datsamparea), mean = 0, sd = sqrt(sig2ehatsamp))
#   ycondgen
# }

EBP_nonsamp_MC = function(BL, Iareas, coefarea, sig2uhatsamp, betahatsamp, datnonsamparea, sig2ehatsamp, coefunit, par_method = "Ori", hlist){
  
  # cond.moments = cond.mu.v(betahatsamp,sig2uhatsamp, sig2ehatsamp, xbaris, ybaris, gammahatis)
  # T.uhatis = cond.moments[[1]] 
  # T.vcondus = cond.moments[[2]]
  
  if(par_method == "Ori"){
    yNs <- replicate(BL, gendata_nonsamp_Inversion(Iareas,coefarea,sig2uhatsamp, betahatsamp,datnonsamparea, sig2ehatsamp,coefunit[2])) #nrow(datsamparea)*BL
    
  }else{
    yNs <- replicate(BL, gendata_nonsamp_Inversion(Iareas,coefarea,sig2uhatsamp, betahatsamp,datnonsamparea, sig2ehatsamp,coefunit[2])) #nrow(datsamparea)*BL
    
  }
  # yNs[datsamparea$Iunits==1]<-datsamparea$yij[datsamparea$Iunits==1]
  gc()  
  
  ## Calculation EBP ##
  nonsamparea = unique(datnonsamparea$area)
  nh = length(hlist)  
  # theta_rep <- sapply(hlist, function(h) {apply(yNs, 2, function(y){ 
  #   tapply((y), datsamparea$area, h)})}) %>% as_tibble() %>% 
  #   setNames(paste0("h", 1:nh)) %>% 
  #   cbind(area = rep(samparea, times = BL))
  # rm(yNs)
  # gc()
  # 
  theta_rep <- sapply(1:nh, function(i) { # changed for gini coefficient
    # if(i != nh){
      apply(yNs,2 , function(y){tapply((y), datnonsamparea$area, hlist[[i]])})
    # }else{
      # apply(yNs,2 , function(y){tapply(exp(y), datnonsamparea$area, hlist[[i]])})}  # muted because of h7 h8 08192021
  }) %>% as_tibble() %>% 
  setNames(paste0("h", 1:nh)) %>% 
  cbind(area = rep(nonsamparea, times = BL))
  
  rm(yNs)
  gc()
  
  estpssampa <- theta_rep %>% group_by(area) %>% summarise_all(mean) ## EBUP ####
  ## observed conditional variance => estimator of MSE leading term
  M1hat <- theta_rep %>% group_by(area) %>% summarise_all(var)     ## M1hat ####
  gc()
  
  return(list (estpssampa = estpssampa
               , M1hat = M1hat))
}



## True MSE ##

Tmse = function(estpssampa,popsampa,nh){
  Tmse.res = sapply(1:nh, function(i){
    colname=paste0("h",i)
    ((estpssampa[,colname])-popsampa[,colname])^2
  })
  colnames(Tmse.res) = paste0("h",1:nh)
  
  Tmse.res  
  
}

##### editted 08032021 for nonsampled area of PS predictor
EBP_nonsamp_PS = function( datsamp, datnonsamparea, Iareas, betahatsamp, b,sig2ehatsamp, uhatis){
  ####PS estimator
  W = tapply(datsamp$probi,datsamp$area,function(.){sample(.,size = 1)}); wi=1/W; names(wi)<-names(Iareas[Iareas==1])
  PSpred.nonarea <- betahatsamp[1] + betahatsamp[2]*tapply(datnonsamparea$x, datnonsamparea$area, mean)+b*sig2ehatsamp+(sum((wi-1)*uhatis)/sum((wi-1)))
  PSpred.nonarea
}





com.vec= function(Iareas,vec){
  Noareas = Iareas[Iareas==0]; Novec = rep(NA,length(Noareas)); names(Novec)=names(Noareas)
  names(vec) = names(Iareas[Iareas==1])
  allareamean = c(vec,Novec)
  allareamean =allareamean[match(names(Iareas),names(allareamean))]
  allareamean
}


