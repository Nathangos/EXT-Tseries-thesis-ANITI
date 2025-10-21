Hillish_stat<-function(Xi,Eta,k){
  Order_Xi<-order(Xi,decreasing = TRUE)
  Eta_star<-Eta[Order_Xi]
  Ni_vector<-sapply(c(1:k),FUN=compute_Nik,
                    k=k,Eta_star=Eta_star)
  vect_k<-c(1:k)
  Hillish_k<-mean(log(k/vect_k)*log(k/Ni_vector))
  return(Hillish_k)
}
Pickands_stat<-function(Xi,Eta,k,q){
  Order_Xi<-order(Xi,decreasing = TRUE)
  Eta_star<-Eta[Order_Xi]
  Sub_eta_sort<-sort(Eta_star[1:k])
  kalt<-(k/2)
  Num_1<-Sub_eta_sort[floor(q*k)]
  Denom_1<-Num_1
  Denom_2<-Sub_eta_sort[floor(q*kalt)]
  k_second<-floor(kalt)
  # Select sub series of smaller size. 
  Sub_eta2_sort<-sort(Eta_star[1:k_second])
  Num_2<-Sub_eta2_sort[floor(q*kalt)]
  Num_<-Num_1-Num_2
  Denom_<-Denom_1-Denom_2
  return(Num_/Denom_)
}
compute_Nik<-function(Eta_star,k,i){
  end<-k
  beg<-i
  Sub_set<-Eta_star[beg:k]
  Value_ref<-Eta_star[i]
  Mean_ki<-sum(as.numeric(Sub_set<=Value_ref))
  return(Mean_ki)
}
Run_diagnostics_Gamma_G<-function(LIMS_Y,dims_elt_text,Vectors_HTAIL,I,J,
                                  y_lims,Vect_k,NAME_Vars,q){
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_min.png"))
  # Dot product --MRV test
  couple_UV<-Vectors_HTAIL[,c(I,J)]
  Nobs<-nrow(couple_UV)
  # Rank transformation. F--> Pareto
  ################
  couple_UV<-apply(X = couple_UV,MARGIN = 2,FUN = function(x){
    Denom<-Nobs+1-rank(x)
    return(Nobs/Denom)
  })
  Min_Risk_functionals_ALL<-apply(X =  couple_UV,MARGIN = 1,
                                  FUN = min)
  GG_min<-Graphics_estimators_gamma(series = Min_Risk_functionals_ALL,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  dev.off()
  #### Max  
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_max.png"))
  Max_Risk_functionals_ALL<-apply(X =  couple_UV,MARGIN = 1,
                                  FUN = max)
  Graphics_estimators_gamma(series = Max_Risk_functionals_ALL,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  dev.off()
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[I],".png"))
  Graphics_estimators_gamma(series = Vectors_HTAIL[,I],
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  
  dev.off()
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[J],".png"))
  obj<-Graphics_estimators_gamma(series = Vectors_HTAIL[,J],
                                 vect_k =Vect_k,
                                 Title_graphic =" ",
                                 dims_elt_text = dims_elt_text,
                                 y_lims = LIMS_Y)
  dev.off()
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_max_",
                    NAME_Vars[I],"_",NAME_Vars[J],".png"))
  Max_Risk_functionals<-apply(X =  couple_UV[,c(I,J)],MARGIN = 1,
                              FUN = max)
  Graphics_estimators_gamma(series = Max_Risk_functionals,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  dev.off()
  
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_min_",
                    NAME_Vars[I],"_",NAME_Vars[J],".png"))
  Min_Risk_functionals<-apply(X =  couple_UV[,c(I,J)],MARGIN = 1,
                              FUN = min)
  Graphics_estimators_gamma(series = Min_Risk_functionals,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  dev.off()
  
  # Confidence bands (ML) --------------------------------------------------------
  ####################
  pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_min_",
                    NAME_Vars[I],"_",NAME_Vars[J],"_conf_band.png"))
  ML_found<-sapply(X =Vect_k ,FUN = function_ML_extRemes,
                   data_d=Min_Risk_functionals,typeML="GP",
                   NB_years=NULL)
  fonction_MLplot_resume(resultatML = ML_found,
                          vecteur_k =Vect_k,nom_variable = " ",
                         lims_Y =LIMS_Y,dims_elt_text=dims_elt_text)
  dev.off()
  
  Rad<-apply(X = couple_UV,MARGIN = 1,FUN = sum)
  Theta<-couple_UV[,1]/Rad
  Hillish_MRV<-sapply(X = Vect_k,FUN = Hillish_stat,
                           Xi=Rad,Eta=Theta)
  Pickandish_MRV<-sapply(X = Vect_k,FUN = Pickands_stat,
                      Xi=Rad,Eta=Theta,q=q)
  png(file = paste0("graphiques_MV/asymp_dependencies/",
                    NAME_Vars[I],"_",NAME_Vars[J],
                    "_Hil_Pick_MRV.png"))
  par(mfrow=c(1,2))
  plot(Vect_k,Hillish_MRV,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value MRV",cex.lab=1.5)
  plot(Vect_k,Pickandish_MRV,type="l",
       xlab="Number of exceedances",
       ylab="Pickandish value MRV",cex.lab=1.5)
  par(mfrow=c(1,1))
  dev.off()
  # Dot product --HRV test ---------------------------------------------------------------------
  ################
  Indexes_<-sapply(X = c(1:nrow(couple_UV)),
                          FUN = function(x){
                            return(ifelse(couple_UV[x,2]>couple_UV[x,1],
                                          yes = 1,
                                          no = 0))
                          })
  
  Indexes_theta_1<-which(Indexes_==1)
  Xi_<-Min_Risk_functionals[Indexes_theta_1]
  Sub_couple<-couple_UV[Indexes_theta_1,]
  Eta_<-Sub_couple[,2]/Sub_couple[,1]
  ## Theta_1
  Hillish_values_1<-sapply(X = Vect_k,FUN = Hillish_stat,
                         Eta=Eta_,Xi=Xi_)
  Hillish_values_2<-sapply(X = Vect_k,FUN = Hillish_stat,
                           Eta=(-1)*Eta_,Xi=Xi_)
  ## Theta_2
  Xi_2<-Min_Risk_functionals[-Indexes_theta_1]
  Sub_couple<-couple_UV[-Indexes_theta_1,]
  Eta_2<-Sub_couple[,1]/Sub_couple[,2]
  Hillish_values_minus_1<-sapply(X = Vect_k,FUN = Hillish_stat,
                           Eta=Eta_2,Xi=Xi_2)
  Hillish_values_minus_2<-sapply(X = Vect_k,FUN = Hillish_stat,
                           Eta=(-1)*Eta_2,Xi=Xi_2)
  png(file = paste0("graphiques_MV/asymp_dependencies/",
                    NAME_Vars[I],"_",NAME_Vars[J],"_Hillish_HRV.png"))
  par(mfrow=c(2,2))
  plot(Vect_k,Hillish_values_1,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5)
  plot(Vect_k,Hillish_values_2,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5)
  plot(Vect_k,Hillish_values_minus_1,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5)
  plot(Vect_k,Hillish_values_minus_2,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5)
  par(mfrow=c(1,1))
  dev.off()
  
  ### Pickands
  png(file = paste0("graphiques_MV/asymp_dependencies/",
                    NAME_Vars[I],"_",NAME_Vars[J],"_Pick_HRV.png"))
  par(mfrow=c(1,2))
  Pickands_1<-sapply(X = Vect_k,FUN=Pickands_stat,
                     q=q,Eta=Eta_,Xi=Xi_)
  Pickands_2<-sapply(X = Vect_k,FUN=Pickands_stat,
                     q=q,Eta=Eta_2,Xi=Xi_2)
  plot(Vect_k,y = Pickands_1,type="l",xlab="Number of exceedances",
       ylab="Pichandish value",cex.lab=1.5)
  plot(Vect_k,y = Pickands_2,type="l",xlab="Number of exceedances",
       ylab="Pichandish value",cex.lab=1.5)
  par(mfrow=c(1,1))
  dev.off()
}
### Simul Mixture Exponential (Resnick)

### (1.13) from Models with Hidden Regular Variation: Generation
### and Detection
#' One_simul_Resnick_Mixture
#'
#' @param Gamma_g : shape parameter (float). 
#' @param p : float (0,1). Probability of the first class.
#' @param vector_theta : vector(float>0). Parameters of exponential
#' distributions. 
#'
#' @return One simulation from the mixture method.
#' @export
#'
#' @examples
One_simul_Resnick_Mixture<-function(Gamma_g,p,vector_theta){
  
  U<-sample(c(1,2),prob=c(p,1-p))
  Vect_return<-rep(1,length(vector_theta))
  Vect_return[U]<-rexp(n = 1,rate = vector_theta)
  simul_pareto<-(1-runif(n =1))^(-Gamma_g)
  return(Vect_return*simul_pareto)
}
Simul_Resnick_Mixture<-function(M,Gamma_g,p,vector_theta){
  M_simuls<-replicate(n = M,expr =One_simul_Resnick_Mixture(Gamma_g=Gamma_g,
                                                  p=p,
                                                  vector_theta=vector_theta))
  return(M_simuls)
}
### Simulation method
Sample_one_mixture_d<-function(prop_bernoulli,d,RiskFunction){
  PASS<-FALSE
  vect_return<-rep(NA,2)
  while(PASS==FALSE){
    Sample_bern<-sample(c(1:d),size = 1,prob = prop_bernoulli)
    Simul_gumbel<--log(-log(runif(1)))
    vect_return<-rep(-Inf,d)
    vect_return[Sample_bern]<-Simul_gumbel
    
    # Rejection_sampling ------------------------------------------------------
    U<-runif(1)
    Ratio_RF<-(exp(RiskFunction(vect_return))
               /sum(exp(vect_return)))
    if(U<Ratio_RF){
      return(vect_return)
    }
  }
}
Mixture_d<-function(prop_bernoulli, nb_simulations,d,RiskFunction){
  Result_mixtures<-replicate(n = nb_simulations,
                             Sample_one_mixture_d(prop_bernoulli = prop_bernoulli,
                                                  d = d,RiskFunction = RiskFunction))
  return(Result_mixtures)
}