######################################################
# Population & Sampling process
######################################################

# This code is modified from functions_MSE.R
# Parameter Estimation :

## Generate fixed population parameters : xij, stata, Nis

genppsfixpop <- function(D){
  
  # Generate area size Nis
  # Nis <- round(1000*(0.5 + runif(D))) %>% setNames(1:D); #Area i sample size
  Nis <- rep(100, D) %>% setNames(1:D)
  # Generate covariate xij
  # tis <- 1 + 3*round(1:D - D/3*round(3/D*(1:D)) )/40 # V
  # names(tis) <- 1:D
  
  areafacpop <- rep(1:D, Nis)
  
  # xij <- tis + eij
  xij <-  runif(sum(Nis), 0, 1)
  
  # Generate stratindpop.
  N_Stratum = round(D/3)
  
  stratindarea <-c( rep(1, N_Stratum),
                    rep(2 , N_Stratum),
                    rep(3, (D - 2*N_Stratum)) ) %>% setNames(1:D)   # names(stratindarea) <- 1:D

  # stratindarea -> stratindpop
  stratindpop <- stratindarea[areafacpop]
  
  list(Nis = Nis, areafacpop = areafacpop, xij = xij, stratindarea = stratindarea, stratindpop = stratindpop)
}

#Step 1-2 

# Generate population yij

genpopPS <- function(D, sig2u, sig2e, beta0, beta1 , Nis, areafacpop, xij,trc=2.5){
  
  # Generate population area random effects : usD
  usD <- rnorm(D, mean = 0, sd = sqrt(sig2u))
  usD[abs(usD) > trc*sqrt(sig2u)] <- sign(usD[abs(usD) > trc*sqrt(sig2u)] )*trc*sqrt(sig2u)
  
  # Generate individual errors : eij
  eij <- rnorm(sum(Nis), mean = 0, sd = sqrt(sig2e))
  eij[abs(eij)> trc*sqrt(sig2e)] <- sign(eij[abs(eij) > trc*sqrt(sig2e)])*trc*sqrt(sig2e)
  
  # Generate response varrables : yij
  yij <- beta0  + beta1*xij + usD[areafacpop] + eij
  
  list(usD = usD, eij = eij, yij = yij)
}

## Generate the sample based on two-step sampling::UPsystematic

gensampPS <- function(D, usD, eij, yij, sig2u, sig2e, stratindarea, stratindpop, Nis,areafacpop){
  # Generate pi
  zi <- round(1000*exp(-usD/8/sqrt(sig2u)))
  ziUs <- tapply(zi, stratindarea, sum)
  
  # eg. D = 50; probi <- 10*zi/(ziUs[stratindarea])
  n_area_Stratum = D*0.6/3
  probi <- n_area_Stratum*zi/(ziUs[stratindarea])
  
  # Generate area indicator VVV
  Iareas <- unlist(lapply(unique(stratindarea), function(i){ sampling::UPsystematic(probi[stratindarea == i])} ))
  names(Iareas) = 1:D
  samparea <- c(1:D)[Iareas == 1]
  
  # Generate pij
  deltaij <- rnorm(sum(Nis), 0, 1)
  zij <- exp( -(usD[areafacpop] + eij)/sqrt(sig2e)/3 + deltaij/15)
  zijUs <-  tapply(zij, areafacpop , sum)
  nis <- rep(5, length(stratindpop))
  nis[stratindpop == 2] <- 10; nis[stratindpop == 3] <- 15
  probij <- nis*zij/(zijUs[areafacpop])
  
  #Generate unit indicator
  Iunits <- unlist(lapply(1:D, function(i){sampling::UPsystematic(probij[areafacpop == i])} )) # VVVVV
  Iunits[!(areafacpop %in% samparea)] <- 0
  
  list(probi = probi, Iareas = Iareas, samparea = samparea,  probij = probij, Iunits = Iunits)
}

# Step 3 - 1 - 2

# Generate yij of sampled area only
simyareasamp <- function(uhatis, vcondus, betahatsamp, datsamp, sig2ehatsamp){
  
  ucondgen <- rnorm(length(uhatis), mean = uhatis, sd = sqrt(vcondus))
  names(ucondgen) <- names(uhatis)
  ycondgen <- betahatsamp[1] + betahatsamp[2]*datsamp$x + ucondgen[as.character(datsamp$area)] + rnorm(nrow(datsamp), mean = 0, sd = sqrt(sig2ehatsamp))
  ycondgen
  
}

######################################################

