Summary_MSE.dfs = function(named_vec){
  # named_vec = avg.vec
  T.name = names(named_vec)
  val = named_vec
  
  MSE_T = grep("MSE",T.name)
  method = gsub("(.*)_h.*","\\1",T.name)
  h_id = gsub(".*_(h.*)", "\\1",T.name)

  long.form = data.frame(h_id=h_id, method= method,val = val) %>% 
    mutate(row_ID = 1:length(named_vec), MSE_T = ifelse(row_ID %in% MSE_T, 1, 0))  %>%
    dplyr::select(MSE_T, h_id, method, val)
  Methods = c("h_id","TMSE","MSE.add","MSE.mult", "MSE.hm","MSE.comp","MSE.noBC", "MSE.ML")#, "MSE.BG1", "MSE.BG2") # changed 08192021
  MSE.df = filter(long.form, MSE_T ==1)[,-1] %>% spread(.,method,val) %>% dplyr::select(one_of(Methods))
  # cbind(h_id = Pred.df[,1], Pred.df[,-1]-Pred.df[,"TR"])
  
  Bias.df = MSE.df
  Bias.df[,-1] = ((Bias.df[,-1]-Bias.df[,"TMSE"])/Bias.df[,"TMSE"])*100
  Bias.df = Bias.df[,-which(colnames(Bias.df)=="TMSE")]
  Bias.df$BEST = colnames(Bias.df)[-1][apply(abs(Bias.df[,-1]),1,which.min)]
  ##########################################################
  M1_T = grep("M1", T.name)
  long.form.M1 = long.form %>% 
                 mutate(row_ID = 1:length(named_vec), M1_T = ifelse(row_ID %in% M1_T, 1, 0))  %>%
                 dplyr::select(M1_T, h_id, method, val)
  M1_method = c("T.M1.bt", "M1.bt","M1hat.bt")
  M1.df = filter(long.form.M1, M1_T ==1)[,-1] %>% spread(.,method,val) %>% 
    dplyr::select(one_of(M1_method))%>% setnames(c("True_M1","M1hat","M1hat_boot"))
  M1.df$True_ratio = M1.df$M1hat / M1.df$True_M1
  M1.df$Est_ratio = M1.df$M1hat_boot / M1.df$M1hat
  M1.df$True_diff = M1.df$M1hat - M1.df$True_M1
  M1.df$Est_diff = M1.df$M1hat_boot - M1.df$M1hat
  
  # # MSE summary h1:
  # summary_h1 = matrix(c(avg.vec[paste0(c("EBP","TMSE","MSE.add","MSE.mult","MSE.hm",
  #                                        "PS","Tmse.PS","mse.PS","mse.PS","mse.PS"),"_h1")]), nrow = 5, ncol = 2, byrow = F) %>% data.frame(.)
  # rownames(summary_h1) = c("predictor","True_MSE","MSE_est.add","MSE_est.mult","MSE_est.hm")
  # colnames(summary_h1) = c("EBP","PS")
  # 
  # ##### confidence Interval

  # CI90_T = grep("90", T.name)
  # long.form.90 = long.form %>%
  #   mutate(row_ID = 1:length(named_vec), CI90_T = ifelse(row_ID %in% CI90_T, 1, 0))  %>%
  #   dplyr::select(CI90_T, h_id, method, val) %>% filter(CI90_T ==1) %>% spread(.,method,val) %>% .[,-1]

  CI95_T = grep("95", T.name)
  long.form.95 = long.form %>%
    mutate(row_ID = 1:length(named_vec), CI95_T = ifelse(row_ID %in% CI95_T, 1, 0))  %>%
    dplyr::select(CI95_T, h_id, method, val) %>% filter(CI95_T ==1) %>% spread(.,method,val) %>% .[,-1]
  # 
  # CI99_T = grep("99", T.name)
  # long.form.99 = long.form %>%
  #   mutate(row_ID = 1:length(named_vec), CI99_T = ifelse(row_ID %in% CI99_T, 1, 0))  %>%
  #   dplyr::select(CI99_T, h_id, method, val) %>% filter(CI99_T ==1) %>% spread(.,method,val) %>% .[,-1]
  # # 
  # 
  #  
  return(list(Bias.df = Bias.df,
              M1.df = M1.df,
              # summary_h1 = summary_h1,
              # CI90 = long.form.90,
              CI95 = long.form.95
              # ,   CI99 = long.form.99
              ))
  
}

