Hillish_stat<-function(Xi,Eta,k){
  Order_Xi<-order(Xi,decreasing = TRUE)
  Eta_star<-Eta[Order_Xi]
  # Ni_vector<-sapply(c(1:k),FUN=compute_Nik,
  #                   k=k,Eta_star=Eta_star)
  Ni_vector<-rank(Eta_star[1:k], ties.method = "average")
  vect_k<-c(1:k)
  Hillish_k <- mean(sapply(1:k, function(j) log(k / Ni_vector[j]) * log(k / j)))
  #Hillish_k<-mean(log(k/vect_k)*log(k/Ni_vector))
  return(Hillish_k)
}
Pickands_stat<-function(Xi,Eta,k,q){
  Order_Xi<-order(Xi,decreasing = TRUE)
  Eta_star<-Eta[Order_Xi]
  Sub_eta<-sort(Eta_star[1:k])
  kalt<-(k/2)
  supkalt<-ceiling(q*kalt)
  Num_1<-Sub_eta[ceiling(q*k)]
  Denom_2<-Sub_eta[supkalt]
  k_second<-ceiling(kalt)
  
  # Select sub series of smaller size. 
  Sub_eta2<-sort(Eta_star[1:k_second])
  Num_2<-Sub_eta2[supkalt]
  Num_<-Num_1-Num_2
  Denom_<-Num_1-Denom_2
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
                                  Vect_k,NAME_Vars,q,convert_HTAIL){
  # Dot product --MRV test
  Nobs<-nrow(Vectors_HTAIL)
  # Rank transformation. F--> Pareto
  ################
  if((convert_HTAIL)==TRUE){
    Vectors_HTAIL<-apply(X =Vectors_HTAIL,MARGIN = 2,FUN = function(x){
      Denom<-Nobs+1-rank(x)
      return(Nobs/Denom)
    })
  }
  if(length(I)==1){
    Name_i<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[I],".png")
    GGi<-Graphics_estimators_gamma(series = Vectors_HTAIL[,I],
                              vect_k =Vect_k,
                              Title_graphic =" ",
                              dims_elt_text = dims_elt_text,
                              y_lims = LIMS_Y)
    ggsave(filename = Name_i,plot =GGi,width=8,
           height=6)
  }
  if(length(J)==1){
    Name_j<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[J],".png")
    obj<-Graphics_estimators_gamma(series = Vectors_HTAIL[,J],
                                   vect_k =Vect_k,
                                   Title_graphic =" ",
                                   dims_elt_text = dims_elt_text,
                                   y_lims = LIMS_Y)
    ggsave(filename = Name_j,plot =obj,width=8,
           height=6)
  }
  L_i<-length(NAME_Vars[I])
  if(L_i>1){
    First_i<-NAME_Vars[I][1]
    for(j in c(2:L_i)){
      First_i<-paste0(First_i,"_",NAME_Vars[I][j])
    }
  }else{
    First_i<-NAME_Vars[I]
  }
  L_j<-length(NAME_Vars[J])
  if(L_j>1){
    Second_j<-NAME_Vars[J][1]
    for(j in c(2:L_j)){
      Second_j<-paste0(Second_j,"_",NAME_Vars[J][j])
    }
  }else{
    Second_j<-NAME_Vars[J]
  }
  Name_max<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_max_",
                   First_i,"_",Second_j,".png")
  Max_Risk_functionals<-apply(X =  Vectors_HTAIL[,c(I,J)],MARGIN = 1,
                              FUN = max)
  GG_max<-Graphics_estimators_gamma(series = Max_Risk_functionals,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  ggsave(filename = Name_max,plot = GG_max,width=8,
         height=6)
  
  Name_min<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_min_",
                   First_i,"_",Second_j,".png")
  Min_Risk_functionals<-apply(X =  Vectors_HTAIL[,c(I,J)],MARGIN = 1,
                              FUN = min)
  GG_min<-Graphics_estimators_gamma(series = Min_Risk_functionals,
                            vect_k =Vect_k,
                            Title_graphic =" ",
                            dims_elt_text = dims_elt_text,
                            y_lims = LIMS_Y)
  ggsave(filename = Name_min,plot = GG_min,width=8,
         height=6)
  
  # Confidence bands (ML) --------------------------------------------------------
  ####################
  # pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_min_",
  #                   First_i,"_",Second_j,"_conf_band.png"))
  # ML_found<-sapply(X =Vect_k ,FUN = function_ML_extRemes,
  #                  data_d=Min_Risk_functionals,typeML="GP",
  #                  NB_years=NULL)
  # fonction_MLplot_resume(resultatML = ML_found,
  #                         vecteur_k =Vect_k,nom_variable = " ",
  #                        lims_Y =LIMS_Y,dims_elt_text=dims_elt_text)
  # dev.off()
  # 
  couple_UV<-Vectors_HTAIL[,c(I,J)]
  d<-ncol(couple_UV)
  if(d==2){
    Rad<-apply(X =couple_UV,MARGIN = 1,FUN = sum)
    Theta<-couple_UV[,1]/Rad
    Hillish_MRV1<-sapply(X = Vect_k,FUN = Hillish_stat,
                        Xi=Rad,Eta=Theta)
    Hillish_MRV2<-sapply(X = Vect_k,FUN = Hillish_stat,
                         Xi=Rad,Eta=-Theta)
    Pickandish_MRV<-sapply(X = Vect_k,FUN = Pickands_stat,
                           Xi=Rad,Eta=Theta,q=q)
    png(file = paste0("graphiques_MV/asymp_dependencies/",
                      NAME_Vars[I],"_",NAME_Vars[J],
                      "_Hil_Pick_MRV.png"))
    par(mfrow=c(1,3))
    plot(Vect_k,Hillish_MRV1,type="l",
         xlab="Number of exceedances",
         ylab="Hillish value MRV",cex.lab=1.5)
    plot(Vect_k,Hillish_MRV2,type="l",
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
  
}
Simple_MRV_HRV_analysis<-function(couple_UV,root_graphics,
                                  Vect_k,l_name,q){
  Min_Risk_functionals<-apply(X = couple_UV,
                              MARGIN = 1,
                              FUN = min)
  ### MRV
  Rad<-apply(X =couple_UV,MARGIN = 1,FUN = sum)
  Theta<-couple_UV[,1]/Rad
  Hillish_MRV1<-sapply(X = Vect_k,FUN = Hillish_stat,
                       Xi=Rad,Eta=Theta)
  Hillish_MRV2<-sapply(X = Vect_k,FUN = Hillish_stat,
                       Xi=Rad,Eta=-Theta)
  Pickandish_MRV<-sapply(X = Vect_k,FUN = Pickands_stat,
                         Xi=Rad,Eta=Theta,q=q)
  png(file = paste0(root_graphics,
                    l_name[1],"_",l_name[2],
                    "_Hil_Pick_MRV.png"))
  par(mfrow=c(1,3))
  plot(Vect_k,Hillish_MRV1,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value MRV",cex.lab=1.5,
       main="Plus case")
  plot(Vect_k,Hillish_MRV2,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value MRV",cex.lab=1.5,
       main="Minus case")
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
  png(file = paste0(root_graphics,
                    l_name[1],"_",l_name[2],"_Hillish_HRV.png"))
  par(mfrow=c(2,2))
  plot(Vect_k,Hillish_values_1,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5,
       main=expression(X[2]>X[1]~" (plus)"))
  plot(Vect_k,Hillish_values_2,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5,
       main=expression(X[2]>X[1]~" (minus)"))
  
  plot(Vect_k,Hillish_values_minus_1,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5,
       main=expression(X[2]<X[1]~" (plus)"))
  plot(Vect_k,Hillish_values_minus_2,type="l",
       xlab="Number of exceedances",
       ylab="Hillish value",cex.lab=1.5,
       main=expression(X[2]<X[1]~" (minus)"))
  par(mfrow=c(1,1))
  dev.off()
  
  ### Pickands
  png(file = paste0(root_graphics,
                    l_name[1],"_",l_name[2],"_Pick_HRV.png"),
      height=600,width=800)
  par(mfrow=c(1,2))
  Pickands_1<-sapply(X = Vect_k,FUN=Pickands_stat,
                     q=q,Eta=Eta_,Xi=Xi_)
  Pickands_2<-sapply(X = Vect_k,FUN=Pickands_stat,
                     q=q,Eta=Eta_2,Xi=Xi_2)
  plot(Vect_k,y = Pickands_1,type="l",xlab="Number of exceedances",
       ylab="Pichandish value",cex.lab=1.5,
       main=expression(X[2]>X[1]))
  plot(Vect_k,y = Pickands_2,type="l",xlab="Number of exceedances",
       ylab="Pichandish value",cex.lab=1.5,
       main=expression(X[2]<X[1]))
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