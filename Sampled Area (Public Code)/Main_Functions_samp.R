rm(list = ls(all = TRUE))

# Sources :
pkgs <- c("lme4","sampling","data.table","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr","reldist")
sapply(pkgs, library, character.only = TRUE, quietly = TRUE)


dir <- "./Cho2024/OptimalPred2024/Sampled Area (Public Code)"
source(file.path(dir, "Summaryfn_MSE_0819.R"))
source(file.path(dir, "fns_population_gen_sampling_0623.R"))
source(file.path(dir, "fns_datasets_0819.R"))
source(file.path(dir, "fns_predictors_0819.R"))
##################################################################################################
s = 1; T.seed = 1;

# Scenarios :
sig2e = (0.3)^2

Scenario = expand.grid( D =  50,
												r.sigma = c(0.5, 1, 2, 3),
												beta0 = 5,
												beta1 = 0.1
) %>% arrange(D, r.sigma, beta0)

D <- Scenario[s,"D"]
sig2u <- sig2e*(Scenario[s,"r.sigma"])^2
beta0 <- Scenario[s,"beta0"]; beta1 <- Scenario[s,"beta1"]
r.sigma <- Scenario[s,"r.sigma"]

#Generate fixed paramters
set.seed(1)
fixpopPS <- genppsfixpop(D)
Nis <- fixpopPS[["Nis"]]; areafacpop <- fixpopPS[["areafacpop"]]; xij <- fixpopPS[["xij"]]; 
stratindarea <- fixpopPS[["stratindarea"]];stratindpop <- fixpopPS[["stratindpop"]]
head(xij)
T.Method = "Ori"
T.variablenames = c("TR", "EBP", "TMSE")
T.par = c("h1","h2","h3","h4","h5","h6") 

variable_names = expand.grid(T.variablenames,T.par) %>% apply(.,1,function(.){paste0(.,collapse = "_")})

##################
## Main Function:
## Return length(variable_names) * D, where only sampled area columns have the returning values:
##################
set.seed(T.seed)
sub.fct.0623(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, T.Method, variable_names,stratindarea, stratindpop)
########################


###############################################################
## Main Function:
### Return length(variable_names) * D matrix, where each column is area.
## only sampled area has values
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
		ind <- ifelse(exp(y) < 155, 1, 0)   ## changed from 90 to 100 : 08162021
		val <- (1-exp(y)/155)*ind
		mean(val)
	}
	h6 = function(y) {gini(exp(y))}
	
	hlist <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4, h5 = h5, h6 = h6)
	
	nh <- length(hlist)
	
	T.par = paste0("h",1:nh)
	names(hlist) = T.par
	
	BL <- 10;  # MC iterations for constructing predictor
	
	
	########################################
	## 1. Data Generation ##
	########################################
	
	# print(1)
	# NDOE: dat.gen() generates the population data, the sample and also fits
	# the sample model to the sampled data. Based on the fitted model, the
	# proposed predictor is computed.
	generated_data = dat.gen(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij,T.Method, stratindarea, stratindpop)
	
	# NDOE: über dieses assign() werden alle Elemente der Ergebnisliste aus generated_data / dat.gen
	# im caller env zugewiesen -> erlaubt verwendung als Funktionsargumente
	for(i in 1:length(generated_data)){ assign(names(generated_data)[i], generated_data[[i]]) }
	gc()
	
	popsampa<- sapply(1:nh, function(i){ 
		# if (i!=nh){      ## modified to transforme exp(yij) for gini coefficient 07302021
		tapply(datsamparea$yij, datsamparea$area, hlist[[i]])}) %>%
		# } else {  tapply(exp(datsamparea$yij), datsamparea$area,hlist[[i]])}}) %>%
		as.data.frame(.) %>% setNames(paste0("h",1:nh)) 
	
	TR = popsampa; # True parameter for sampled areas (dim: sampled areas * considered parameter) 
	
	################################################
	## 2. EBP calculation of the proposed predictor:
	################################################
  EBP_M1_MC = EBP_MC(BL, betahatsamp, sig2uhatsamp, sig2ehatsamp, xbaris, ybaris, gammahatis, datsamparea, coefunit, par_method = T.Method, hlist)
	EBP = EBP_M1_MC$estpssampa %>% as.matrix(.) %>% .[,T.par]
	
	TMSE = Tmse(EBP,TR,nh) 
	
	
	########################################
	# Save results to variable_names
	########################################
	T_method = gsub("(.*)_h.*","\\1",variable_names) %>% unique(.)
	T_par = paste0("h",1:6)
	
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

##################
## Main Function:
## Return length(variable_names) * D, where only sampled area columns have the returning values:
##################
set.seed(T.seed)
sub.fct.0623(D, sig2u, sig2e, beta0, beta1, Nis, areafacpop, xij, T.Method, variable_names,stratindarea, stratindpop)
########################



