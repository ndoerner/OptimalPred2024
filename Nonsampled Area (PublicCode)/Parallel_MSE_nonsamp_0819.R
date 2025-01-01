rm(list = ls(all = TRUE))

# Sources :
all(sapply(c("lme4","sampling","data.table","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr","reldist"), require, character.only = TRUE))
source("Parallel_MSE_nonsamp_0819.R")
source("fns_population_gen_sampling_0623.R")
source("fns_datasets_nonsamp_0819.R")
source("fns_predictors_nonsamp_0819.R")

#Scenario Parameter configuration:
s = 1; T.seed = 1;

sig2e = (0.3)^2

Scenario = expand.grid( D =  150,
                        r.sigma = c(0.5, 1, 2, 3),
                        beta0 = 5,
                        beta1 = 0.1
) %>% arrange(D, r.sigma, beta0)

# Parameters Generation from the Scenario matrix :
D <- Scenario[s,"D"]
sig2u <- sig2e*(Scenario[s,"r.sigma"])^2
beta0 <- Scenario[s,"beta0"]; beta1 <- Scenario[s,"beta1"]
r.sigma <- Scenario[s,"r.sigma"]

#Generate fixed paramters
set.seed(1)
fixpopPS <- genppsfixpop(D)
Nis <- fixpopPS[["Nis"]]; areafacpop <- fixpopPS[["areafacpop"]]; xij <- fixpopPS[["xij"]]; 
stratindarea <- fixpopPS[["stratindarea"]];stratindpop <- fixpopPS[["stratindpop"]]
T.Method = "Ori"
head(xij)

T.variablenames = c("TR", "EBP", "TMSE")
T.par = c("h1","h2","h3","h4","h5","h6") # updated h8 gini(abs(y))

variable_names = expand.grid(T.variablenames,T.par) %>% apply(.,1,function(.){paste0(.,collapse = "_")})

#######################################################################################################################
set.seed(T.seed)
MC_one_iter = sub.fct.0623(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, T.Method, variable_names,stratindarea, stratindpop)

#######################################################################################################################


###############################################################
## Main Function:
### Return length(variable_names) * D matrix, where each column is area.
###############################################################


sub.fct.0623 = function(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, T.Method, variable_names,stratindarea, stratindpop){
  
  ########################################
  ## 0. Parameter Definition:
  ########################################
  h1 = function(y) {mean(y)}
  h2 = function(y) {mean(exp(y))}
  h3 = function(y) {quantile(y, probs = 0.25)} #Q1
  h4 = function(y) {quantile(y, probs = 0.75)} #Q3
  h5 = function(y) { #PI
    ind <- ifelse(exp(y) < 155, 1, 0)  
    val <- (1-exp(y)/155)*ind
    mean(val)
  }
  h6 = function(y) {gini(exp(y))}
  
  hlist <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4, h5 = h5, h6 = h6)
  
  nh <- length(hlist)
  
  T.par = paste0("h",1:nh)
  names(hlist) = T.par

  BL <- 250; # MC iterations for constructing predictor
  
  ########################################
  ## 1. Data Generation ##
  ########################################
  # print(1)
  generated_data = dat.gen(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij,T.Method, stratindarea, stratindpop)
  for(i in 1:length(generated_data)){ assign(names(generated_data)[i], generated_data[[i]]) }
  gc()

  ## true value ####
  popsampa<- sapply(1:nh, function(i){ 
    tapply(datnonsamparea$yij, datnonsamparea$area, hlist[[i]])})%>%
    as.data.frame(.) %>% setNames(paste0("h",1:nh)) 
  TR = popsampa;
  
  ################################################
  ## 2. EBP calculation of the proposed predictor:
  ################################################
  ## Generate yij where j is ri ###
  EBP_M1_MC = EBP_nonsamp_MC_CI(BL, Iareas, coefarea, sig2uhatsamp, betahatsamp, datnonsamparea, sig2ehatsamp, coefunit, par_method = "Ori", hlist, popsampa)
  EBP = EBP_M1_MC$estpssampa %>% as.matrix(.) %>% .[,T.par] 
  
  TMSE = Tmse(EBP,TR,nh) # Calculate True MSE
  
  ########################################
  # Save results to variable_names
  ########################################
  
  T_method = gsub("(.*)_h.*","\\1",variable_names)%>% unique(.)
  T_par = paste0("h",1:nh)
  
  for(i in T_method){
    for(j in T_par){
      Text = paste0(i,"_",j," = ",i,"[,'",j,"']")
      eval(parse(text = Text))
    }
  }
  
  res_df = lapply(variable_names, function(.){eval(parse(text = .))}) %>% Reduce("rbind",.)
  null_df = matrix(NA, nrow = nrow(res_df), ncol = D)
  colnames(null_df) = as.character(1:D)
  null_df[,rownames(TR)] = res_df
  
  return(null_df)
}  
 
  
  
