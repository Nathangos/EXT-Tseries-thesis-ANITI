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
Run_diagnostics_Gamma_G<-function(dims_elt_text,Vectors_HTAIL,
                                  Vect_k,q,root_export,YLIM_MV){
  Mat_combs<-t(utils::combn(x = c(1:ncol(Vectors_HTAIL)),
                        m = 2))
  list_gg<-list()
  Names<-colnames(Vectors_HTAIL)
  for(j in c(1:nrow(Mat_combs))){
    pair<-Mat_combs[j,]
    pair_names<-Names[pair]
    Sub_vector<-Vectors_HTAIL[,pair]
    Minj<-apply(X = Sub_vector,MARGIN = 1,FUN = min)
    GG<-Graphics_estimators_gamma(series = Minj,
                              vect_k = Vect_k,
                              Title_graphic = " ",
                              dims_elt_text = dims_elt_text,
                              SELECT_estim = c("ML_gamma"))
    DF_j<-GG$data
    DF_j$pair_vars<-rep(paste0("Min(",pair_names[1],",",
                               pair_names[2],")"),nrow(DF_j))
    list_gg[[j]]<-DF_j
  }
  Df_all<-do.call(rbind.data.frame,list_gg)
  COLS_chosen<-c("gam_ref"="red",
                 "ML_gamma"="blue",
                 "confidence_band"="darkblue",
                 "threshold"="red")
  LINETYPE_chosen<-c("gam_ref"=3,
                     "Hill_gamma"=2,
                     "ML_gamma"=4,
                     "confidence_band"=5,
                     "threshold"=6)
  Df_all$pair_vars<-sapply(Df_all$pair_vars,
                           FUN =function(x){
                             Fct_correct_name(x = x,target = "Surcote",
                                              replacement = "Surge")})
  GG_shape_AD<-ggplot2::ggplot(data=Df_all,aes(x=number_excesses,y =gamma_estimed,color=source,
                                      group=interaction(source),
                                      linetype=source))+
    geom_line()+
    facet_wrap(~pair_vars)+
    ylab(expression(gamma))+
    xlab("Number of exceedances")+
    ylim(YLIM_MV)+
    geom_ribbon(mapping = aes(ymin=bound_inf,ymax=bound_sup,
                              col="confidence_band",
                              linetype = "confidence_band"),alpha=0.15,
                fill="grey")+
    geom_hline(aes(yintercept=Shape_chosen,col="gam_ref",
                   linetype="gam_ref"))+
    labs(col="Legend", linetype = "Legend")
  LABELS_shape<-names(COLS_chosen)
  LEGEND_shape<-sapply(LABELS_shape,
                       function(x){
                         if(x=="gam_ref"){
                           return(latex2exp::TeX("$\\gamma_{0}$"))
                         }else{
                           return(x)
                         }
                       })
  GG_shape_new<-GG_shape_AD+
    scale_linetype_manual(values = LINETYPE_chosen,
                          labels=LEGEND_shape)+
    scale_color_manual(values = COLS_chosen,
                       labels=LEGEND_shape)+
    theme_bw()+
    theme(axis.title=element_text(size=dims_elt_text[1]),
          legend.title = element_text(size=dims_elt_text[2]),
          legend.text=element_text(size=dims_elt_text[3]),
          axis.text=element_text(size=dims_elt_text[4]))+
    theme(legend.position="bottom",
          legend.direction = "horizontal")
  ggsave(filename = paste0(root_export,"/ADiag_shape_d=",ncol(Vectors_HTAIL),
                           ".png"),
         plot=GG_shape_new,width = 8,height = 6)
  return(NA)
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
Window_std_per_threshold<-function(k_end,bandwidth_h,series_orig){
  Inds_taken<-k_end-bandwidth_h
  ### take results found with higher threshold--> go backwards
  sub_series<-series_orig[Inds_taken:k_end]
  return(sd(sub_series))
}
Stability_criterion_choice_threshold<-function(series_sorted_stats,k_max,
                                               bandwidth_h){
  ind_beg<-bandwidth_h
  Inds_candidates<-c(ind_beg:k_max)
  Vect_std_window<-sapply(Inds_candidates,
         FUN = Window_std_per_threshold,
         series_orig=series_sorted_stats,
         bandwidth_h=bandwidth_h)
  ### Detect local minimal for the std vector
  Differences<-diff(Vect_std_window)
  L_end<-length(Differences)-1
  ### The previous index gives a higher value
  Index_1<-which(Differences[1:L_end]<0)+1
  ### The next index gives a higher value
  d2<-Differences
  Index_2<-which(d2>=0)
  ind_min_local<-intersect(Index_1,Index_2)
  
  ### Determine the local minimum giving a value smaller
  ### than the mean
  ref_val<-mean(Vect_std_window)
  All_candidates<-Vect_std_window[ind_min_local]
  sub_ind_candidates<-which(All_candidates<ref_val)
  index_special<-ind_min_local[sub_ind_candidates]
  
  ### Since we go back to the statistic--> we need to
  ### recover the original index--> add bandwidth_h
  Anchor<-index_special[1]+bandwidth_h
  end_anchor<-index_special[1]
  Inds_sub<-c(Anchor:end_anchor)
  ### select the observations close to the chosen beta
  Sub_stat<-series_sorted_stats[Inds_sub]
  
  ### pick among these observations the one closer to the median
  Med_sub_stat<-median(Sub_stat)
  Ind_beta_star<-which.min(abs(Sub_stat-Med_sub_stat))
  beta_star<-Inds_sub[Ind_beta_star]
  ### (graphics) recover the associated standard deviation (if available)
  newbeta_star<-ifelse(beta_star>bandwidth_h,
                    yes = beta_star-bandwidth_h,
         no = NA)
  print(c(beta_star,newbeta_star))
  STAT_beta<-ifelse(beta_star>bandwidth_h,
                    yes=Vect_std_window[newbeta_star],
                    no=NA)
  return(list("evol_std"=Vect_std_window,
              "beta_found"=newbeta_star,
              "candidates"=index_special,
              "stat_candidates"=Vect_std_window[index_special],
              "stat_beta"=STAT_beta))
}
Concatenate_entries_same_key<-function(list_,key_target){
  vector_results<-c(sapply(which(names(list_)==key_target),
           FUN=function(x){
             return(list_[[x]])
           }))
  return(vector_results)
}
