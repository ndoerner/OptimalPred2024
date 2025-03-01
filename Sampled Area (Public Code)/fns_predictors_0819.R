# This code is modified from functions_MSE.R
# This code is modified from fns_predicots_0623.R : parameters are changed!!!! ex. Gini coefficient
# This code is modified from fns_predicots_0715.R : use the  "exp(yij)" for Gini coefficient
# This code is modified from fns_predicots_0805.R : changed in EBP_MC for including function "h8" : gini(abs(y)) 08192021


## EBP calculation by Methods:

# NDOE: Erwartungswert und Varianz fuer NV in MC des EBP
cond.mu.v = function(betahat, sig2uhatsamp, sig2ehatsamp, xbaris.w, ybaris.w, gammahatis.w ){
  uhatis.w <- gammahatis.w*(ybaris.w - betahat[1] - betahat[2]*xbaris.w)
  vcondus.w <- (1-gammahatis.w)*sig2uhatsamp
  return(list(uhatis.w = uhatis.w
              , vcondus.w = vcondus.w))
  
}

# NDOE: Zuege aus NV/Erzeugung der Bootstrap-Population
gendata <- function(uhatis,vcondus,betahatsamp,datsamparea,sig2ehatsamp,b){
  ucondgen <- rnorm(length(uhatis), mean = uhatis, sd = sqrt(vcondus))
  names(ucondgen) <- names(uhatis)
  ycondgen <- b*sig2ehatsamp + betahatsamp[1] + betahatsamp[2]*datsamparea$x + ucondgen[as.character(datsamparea$area)] + rnorm(nrow(datsamparea), mean = 0, sd = sqrt(sig2ehatsamp))
  ycondgen
}
# NDOE: wird eigentlich nur genutzt, um den Simulationsteil des EBP durchzuefuehren
EBP_MC = function(BL, betahatsamp, sig2uhatsamp, sig2ehatsamp, xbaris, ybaris, gammahatis, datsamparea, coefunit, par_method = "Ori", hlist){
  
  cond.moments = cond.mu.v(betahatsamp,sig2uhatsamp, sig2ehatsamp, xbaris, ybaris, gammahatis)
  T.uhatis = cond.moments[[1]] 
  T.vcondus = cond.moments[[2]]
  
  if(par_method == "Ori"){
    yNs <- replicate(BL, gendata(T.uhatis , T.vcondus, betahatsamp, datsamparea , sig2ehatsamp, b=coefunit[2])) #nrow(datsamparea)*BL
    
  }else{
    yNs <- replicate(BL, gendata(T.uhatis , T.vcondus, betahatsamp, datsamparea , sig2ehatsamp, b=0)) #nrow(datsamparea)*BL
    
  }
  yNs[datsamparea$Iunits==1]<-datsamparea$yij[datsamparea$Iunits==1]
  gc()  
  
  ## Calculation EBP ##
  samparea = unique(datsamparea$area)
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
      apply(yNs,2 , function(y){tapply((y), datsamparea$area, hlist[[i]])}) 
    # }else{ 
      # apply(yNs,2 , function(y){tapply(exp(y), datsamparea$area, hlist[[i]])})}   # muted to include function "h8" gini(abs(y))
  }) %>% as_tibble() %>%
  setNames(paste0("h", 1:nh)) %>% 
  cbind(area = rep(samparea, times = BL))
  
  rm(yNs)
  gc()
  
  estpssampa <- theta_rep %>% group_by(area) %>% summarise_all(mean) ## EBUP ####
  ## observed conditional variance => estimator of MSE leading term
  M1hat <- theta_rep %>% group_by(area) %>% summarise_all(var)     ## M1hat ####
  rm(theta_rep) # editted 08172021
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

EBP_PS = function(nis, Nis, Iareas, uhatis, betahatsamp, datsamparea, ybaris,xbaris,coefunit,sig2ehatsamp){
  Xi = tapply(datsamparea$x, datsamparea$area, mean);xi = xbaris
  niNi= nis/Nis[Iareas==1]
  PSpred.area <- (1-niNi)*(uhatis+betahatsamp[1]+Xi*betahatsamp[2]+coefunit[2]*sig2ehatsamp)+niNi*(ybaris+betahatsamp[2]*(Xi-xi))
  rm(datsamparea);gc();
  PSpred.area
}

com.vec= function(Iareas,vec){
  Noareas = Iareas[Iareas==0]; Novec = rep(NA,length(Noareas)); names(Novec)=names(Noareas)
  names(vec) = names(Iareas[Iareas==1])
  allareamean = c(vec,Novec)
  allareamean =allareamean[match(names(Iareas),names(allareamean))]
  allareamean
}


