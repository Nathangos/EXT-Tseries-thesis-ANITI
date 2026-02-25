#' Title
#'
#' @param gam_t 
#' @param sig_t 
#' @param Kappa_t 
#' @param p 
#'
#' @return gradient for rl level confidence.Delta method
#' @export
#'
#' @examples
fprime_Rl_extgp<-function(gam_t,sig_t,Kappa_t,p){
  
  ### beta0=shape_param, exp(beta1)=sigma, exp(beta2)=Kappa
  shape_param<-gam_t
  beta_1<-log(sig_t)
  beta_2<-log(Kappa_t)
  rlevel<-as.numeric(mev::qextgp(kappa = Kappa_t,
                                 sigma =  sig_t,
                                 xi = gam_t,type = 1,
                                 p=p))
  ## Derivative|shape
  #u'.v
  grad_1num1<-shape_param*(-log(1-p*exp(-beta_2))*(1-p*exp(-beta_2))^(-shape_param))*exp(beta_1)
  #u.v'
  grad_1num2<-rlevel*shape_param
  grad_1num<-grad_1num1-grad_1num2
  grad_1denom<-shape_param^(2)
  grad_1<-grad_1num/grad_1denom
  
  ## Derivative|scale
  grad_2<-rlevel
  
  ## Derivative|kappa
  Deriv_kappa<-(log(p)*exp(-beta_2))*(p^(exp(-beta_2)))
  grad_3num<-(-shape_param)*(1-p^(exp(-beta_2)))^(-shape_param-1)*Deriv_kappa
  grad_3denom<-shape_param*exp(-beta_1)
  grad_3<-grad_3num/grad_3denom
  
  Vect_grad<-c(grad_1,grad_2,grad_3)
  return(as.numeric(Vect_grad))
}
Fct_correlations<-function(method_corr,df1,df2){
  vect_r<-sapply(X = c(1:ncol(df1)),FUN = function(j,x_1,x_2){
    series_1<-x_1[,j]
    series_2<-x_2[,j]
    return(cor(x = series_1,y = series_2,method = method_corr))
  },x_1=df1,x_2=df2)
  return(vect_r)
}

Resamples_correlations<-function(l_extremes_indus,B,method_corr,l_names){
  N<-nrow(l_extremes_indus[[1]])
  Indexes_B<-sample(x = c(1:N),size = B,replace = TRUE)
  DF1<-l_extremes_indus[[l_names[1]]][Indexes_B,]
  DF2<-l_extremes_indus[[l_names[2]]][Indexes_B,]
  Correlations_sample<-Fct_correlations(method_corr = method_corr ,
                                        df1 =DF1,
                                        df2=DF2)
  
  return(Correlations_sample)
}

######## Consistency measurement

Conf_interval_correlations_2V<-function(df1,df2,niv_conf,method_corr,NB_boot){
  Realisations_values<-replicate(n = NB_boot,Sample_level_correlation(df1=df1,
                                                                      df2=df2,
                                                                      method_corr=method_corr))
  Bound_sup<-apply(X = Realisations_values,MARGIN = 1,FUN = function(x){
    return(quantile(x,(1-(niv_conf/2))))})
  Mean<-apply(X=Realisations_values,MARGIN = 1,FUN = mean)
  Bound_inf<-apply(X = Realisations_values,MARGIN = 1,FUN = function(x){
    return(quantile(x,niv_conf/2))})
  return(cbind(Bound_inf,Mean,Bound_sup))
}
Conf_interval_ext_correlations_2V<-function(df1,df2,niv_conf,method_corr,NB_boot){
  Realisations_values<-replicate(n = NB_boot,Sample_level_correlation(df1=df1,
                                                                      df2=df2,
                                                                      method_corr=method_corr))
  Bound_sup<-apply(X = Realisations_values,MARGIN = 1,FUN = function(x){
    return(quantile(x,(1-(niv_conf/2))))})
  Mean<-apply(X=Realisations_values,MARGIN = 1,FUN = mean)
  Bound_inf<-apply(X = Realisations_values,MARGIN = 1,FUN = function(x){
    return(quantile(x,niv_conf/2))})
  return(cbind(Bound_inf,Mean,Bound_sup))
}
Sample_level_correlation<-function(df1,df2,method_corr){
  Sample_inds<-sample(c(1:nrow(df1)),size =nrow(df1),replace = TRUE )
  Values_obtained<-Fct_correlations(method_corr = "kendall",
                                    df1 =df1[Sample_inds,],
                                    df2=df2[Sample_inds,])
  return(Values_obtained)
}
Cross_Extremo_MV_all<-function(list_simul,list_reality,
                               q_per_time,B,
                               list_names_VAR){
  ### i,j equals
  d<-length(list_names_VAR)
  Full_indexes<- expand.grid(
    row = seq_len(d),
    col = seq_len(d)
  )
  All_C_extremo<-apply(X = Full_indexes,FUN = function(couple_used){
    I<-couple_used[1]
    J<-couple_used[2]
    name_i<-list_names_VAR[I]
    name_gi<-name_i
    if(name_i=="Surcote"){
      name_gi<-"Surge"
    }
    name_j<-list_names_VAR[J]
    name_gj<-name_j
    if(name_j=="Surcote"){
      name_gj<-"Surge"
    }
    
    Name_pair<- paste0("(",name_gi,","
                       ,name_gj,")")
    expr<-expression(gamma[Name_pair]~"value")
    ### do.call etc... to replace in an expr
    ### the variable names by their actual value
    
    Ylab_value<-do.call("substitute", list(expr[[1]], 
                                           list(Name_pair=Name_pair)))
    
    df_cross_extremo_ij<-Cross_Extremo_MV_ij(list_simul = list_simul,
                                             list_reality = list_reality,
                                             name_cond = name_i,name_other=name_j,
                                             q_per_time=q_per_time,B=B)
    df_cross_extremo_ij$category<-rep(Name_pair,
                                      nrow(df_cross_extremo_ij))
    return(df_cross_extremo_ij)
    # GG_gam_ij<-ggplot(data=df_cross_extremo_ij,
    #        aes(x=time_d,y=val_data,col="data"))+
    #   geom_point()+
    #   geom_line()+
    #   geom_point(aes(x=time_d,y=val_sim,
    #                  col="simulations"),pch=2)+
    #   geom_line(aes(x=time_d,y=val_sim,col="simulations"))+
    #   geom_ribbon(mapping = aes(ymin=Qinf,ymax=Qsup,col="confidence_band"),
    #               alpha=0.15,
    #               fill="grey", linetype = "dashed")+
    #   scale_color_manual(values = cols_)+
    #   ylab(Ylab_value)+xlab("Lag h")+
    #   labs(col="Legend")+
    #   theme(axis.title=element_text(size=20),
    #         legend.text=element_text(size=12),
    #         legend.title = element_text(size=13),
    #         axis.text = element_text(size=14))+
    #   scale_linetype_manual("Legend",
    #                         values=c("confidence_band"=2,
    #                                  "mean"=5))
    # return(GG_gam_ij)
  },MARGIN = 1)
  ### Obtained_graphics
  return(All_C_extremo)
}
Cross_Extremo_MV_ij<-function(list_simul,list_reality,name_cond,name_other,
                              q_per_time,B){
  Variable_simul_cond<-list_simul[[name_cond]]
  Variable_reality_cond<-list_reality[[name_cond]]
  ##
  Variable_simul_other<-list_simul[[name_other]]
  Variable_reality_other<-list_reality[[name_other]]
  z<-1
  vect_distances<-c()
  Matrix_pairs<-list()
  L<-ncol(Variable_simul_other)
  for(j in 1:L){
    valeurs_t_plus_h<-j:L
    for (i in valeurs_t_plus_h){
      Matrix_pairs[[z]]<-c(i,j)
      vect_distances<-c(vect_distances,abs(i-j))
      z<-z+1
    }
  }
  Tau_cond<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_reality_cond)
  Tau_other<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_reality_other)
  list_Tau<-list("other"=Tau_other,"cond"=Tau_cond)
  ##
  CROSS_extremo_reality<-cross_extremogram(Matrix_pairs = Matrix_pairs,
                                           var_cond = Variable_reality_cond,
                                           var_other = Variable_reality_other,
                                           l_Tau = list_Tau)
  ### Attention ! We use then the quantile from simulations
  Tau_cond_sim<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_simul_cond)
  Tau_other_sim<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_simul_other)
  list_Tau_sim<-list("other"=Tau_other_sim,
                     "cond"=Tau_cond_sim)
  CROSS_extremo_sim<-cross_extremogram(Matrix_pairs = Matrix_pairs,
                                       var_cond = Variable_simul_cond,
                                       var_other = Variable_simul_other,
                                       l_Tau = list_Tau_sim)
  
  # Bootstrap confidence_intervals ------------------------------------------
  Boot_extremo_results<-replicate(n = B,Resampling_cross_extremo(vect_distances = vect_distances,
                                                                 Matrix_pairs = Matrix_pairs,
                                                                 var_cond = Variable_reality_cond,
                                                                 var_other = Variable_reality_other,
                                                                 l_Tau = list_Tau))
  Extremo_inf<-apply(X =Boot_extremo_results,MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.025))})
  Extremo_sup<-apply(X =Boot_extremo_results,MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.975))})
  df_extremo<-cbind.data.frame(vect_distances,CROSS_extremo_reality,
                               CROSS_extremo_sim)
  colnames(df_extremo)<-c("time_d","reality_extremo","sim_extremo")
  result_delta<-df_extremo %>% 
    group_by(time_d) %>% 
    summarise(val_data=mean(reality_extremo),
              val_sim=mean(sim_extremo))
  result_delta$Qinf<-Extremo_inf
  result_delta$Qsup<-Extremo_sup
  return(result_delta)
}

Resampling_cross_extremo<-function(vect_distances,
                                   Matrix_pairs ,
                                   var_cond,
                                   var_other,
                                   l_Tau){
  Inds_taken<-nrow(var_cond)
  Resample_indexes<-sample(1:Inds_taken,size = Inds_taken,
                           replace = TRUE)
  var_cond_alt<-var_cond[Resample_indexes,]
  var_other_alt<-var_other[Resample_indexes,]
  Cross_boot_extremo<-cross_extremogram(Matrix_pairs =Matrix_pairs,
                                        var_cond = var_cond_alt,
                                        var_other = var_other_alt,
                                        l_Tau = l_Tau)
  df_extremo<-cbind.data.frame(vect_distances,Cross_boot_extremo)
  colnames(df_extremo)<-c("time_d","boot_extremo")
  result_delta_boot<-df_extremo %>% group_by(time_d) %>% summarise(val_boot=mean(boot_extremo))
  return(as.numeric(result_delta_boot$val_boot))
}
F_bounds_per_column<-function(x,Alpha_q){
  return(c(as.numeric(
    quantile(x,
             c(Alpha_q/2,1-(Alpha_q/2)
             )))))
}
cross_extremogram<-function(Matrix_pairs,var_cond,var_other,l_Tau){
  PIST_empirique<-sapply(Matrix_pairs,FUN = cross_pi_s_t,
                         var_cond=var_cond,var_other=var_other,l_Tau=l_Tau)
  return(PIST_empirique)
}
cross_pi_s_t<-function(couple_s_t,var_cond,var_other,l_Tau){
  t<-couple_s_t[1]
  s<-couple_s_t[2]
  v_cond<-l_Tau[["cond"]]
  v_other<-l_Tau[["other"]]
  if(is.null(nrow(var_cond))==FALSE){
    numerateur<-sum(as.numeric((var_cond[,s]>v_cond[s])&
                                 (var_other[,t]>v_other[t])))
    denominateur<-sum(as.numeric(var_cond[,s]>v_cond[s]))
  }
  else{
    numerateur<-sum(as.numeric((var_cond[s]>v_cond[s])&
                                 (var_other[t]>v_other[t])))
    denominateur<-sum(as.numeric(var_cond[s]>v_cond[s]))
  }
  return(numerateur/denominateur)
  
}
F_automatic_ggplot_marg<-function(vect_time_choice, 
                                  l_unit,risk_function,
                                  ind_var_cond,levels_used,
                                  cols_,L_simul,
                                  L_obs,L_nameV,L_whole_obs,
                                  Path_TREND,Path_EXT,hearts){
  obj<-Marg_2d_simul_vs_obs(list_simul = L_simul,
                            list_obs = L_obs,
                            vect_times = vect_time_choice,
                            l_name_2V = L_nameV ,
                            cols_ggplot=cols_,
                            l_unit = l_unit)
  ggsave(filename = paste0(Path_TREND
                           ,"/",Prefix,"_distrib_values_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,
                           opt_used,".png"),
         width=8,height = 6,
         plot = obj)
  
  ## Bivariate return level
  ##### on the original data observations.
  BivRL<-Bivariate_RLevel_simul_vs_obs(list_simul = L_simul,
                                       list_obs = L_whole_obs,
                                       vect_times = vect_time_choice,
                                       l_name_2V =L_nameV ,
                                       cols_ggplot=cols_,
                                       l_unit = l_unit,
                                       levels_used=levels_used,
                                       cols_ = cols_,hearts=hearts)
  ggsave(filename = paste0(Path_EXT
                           ,"/",Prefix,"_BivReturnLevel_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,
                           opt_used,".png"),
         width=8,height = 6,
         plot = BivRL)
  
  ## Cond_ext
  obj_cond<-Marg_2d_cond_ext(list_simul = L_simul,
                             list_obs = L_obs,
                             vect_times = vect_time_choice,
                             l_name_2V =L_nameV ,
                             cols_ggplot=cols_,
                             l_unit = l_unit,
                             risk_function =risk_function,
                             ind_var_cond = ind_var_cond)
  ggsave(filename = paste0(Path_TREND,"/",Prefix,
                           "_cond_distrib_values_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,
                           opt_used,".png"),
         width=8,height = 6,
         plot = obj_cond)
}
#' Title
#'
#' @param gamma 
#' @param sigma_EGPD 
#' @param kappa 
#' @param order_quantiles 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
Ratio_log_model1<-function(gamma,sigma_EGPD,order_quantiles,data){
  
  Xfound<-as.numeric(quantile(data,order_quantiles))
  logemp_fonction<-log(1-(1+gamma*(Xfound/sigma_EGPD))^(-1/gamma))/log(order_quantiles)
  return(logemp_fonction)
}
RL_EGPD_model1<-function(gamma,sigma_EGPD,order_quantiles,Kappa){
  return(mev::qextgp(p = order_quantiles,kappa =Kappa,
                     sigma = sigma_EGPD,xi = gamma,
                     type = 1))
}
#' Marg_2d_simul_vs_obs
#'
#' @param list_simul: list [dataframe]. For each key, a dataframe 
#' of simulated time series. 
#' @param list_obs: list [dataframe].For each key, a dataframe 
#' of simulated time series. 
#' @param vect_times: vector[int]. Vector of analysed times.
#' @param l_name_2V: list[str]. Analysed variables.
#'@param l_unit: list[str]. Unit_used
#'
#' @return GGplot.
#' @export
#'
#' @examples
Marg_2d_simul_vs_obs<-function(list_simul,list_obs,vect_times,l_name_2V,
                               cols_ggplot,l_unit){
  
  Z<-1
  NSimul_for_compar<-nrow(list_simul[[1]])
  NObs_for_compar<-nrow(list_obs[[1]])
  Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
                   ncol =length(vect_times) )
  Obs_ij<-matrix(NA,nrow = NObs_for_compar,
                 ncol =length(vect_times) )
  for(j in c(1:length(list_simul))){
    if(names(list_simul)[j]%in%l_name_2V){
      Simul_namev<-list_simul[[j]][,vect_times[Z]]
      Simul_ij[,Z]<-Simul_namev
      Obs_namev<-list_obs[[j]][,vect_times[Z]]
      Obs_ij[,Z]<-Obs_namev
      Z<-Z+1
    }
  }
  Obs_ij<-as.data.frame(Obs_ij)
  colnames(Obs_ij)<-sapply(c(1:length(l_name_2V)),function(x){
    return(paste0("V_",x))
  })
  Simul_ij<-as.data.frame(Simul_ij)
  colnames(Simul_ij)<-colnames(Obs_ij)
  df_combo<-rbind.data.frame(Obs_ij,Simul_ij)
  df_combo$type<-c(rep("data",nrow(Obs_ij)),
                   rep("simulations",nrow(Simul_ij)))
  exp_type<-expression(x[M]^t~v)
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[1],
                                              x=l_name_2V[1],
                                              v=l_unit[[l_name_2V[1]]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[2],
                                              x=l_name_2V[2],
                                              v=l_unit[[l_name_2V[2]]])))
  GG0<-ggplot(data=Obs_ij,aes(x=V_1,y=V_2,col="data"))+
    geom_point(size=2,pch=19)+
    geom_point(data=Simul_ij,aes(x=V_1,y=V_2,col="simulations"),
               size=0.75,
               pch=17)+
    geom_xsidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    labs(col="Legend")+
    scale_color_manual(values=cols_ggplot)+
    scale_fill_manual(values= cols_ggplot)+
    guides(fill="none")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    xlab(xlab_title)+ylab(ylab_title)
  return(GG0)
  
}

Bivariate_RLevel_simul_vs_obs<-function(list_simul,list_obs,
                                        vect_times,l_name_2V,
                                        cols_ggplot,l_unit,
                                        levels_used,cols_,
                                        hearts){
  Z<-1
  NSimul_for_compar<-nrow(list_simul[[1]])
  NObs_for_compar<-nrow(list_obs[[1]])
  Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
                   ncol =length(vect_times) )
  Obs_ij<-matrix(NA,nrow = NObs_for_compar,
                 ncol =length(vect_times) )
  for(j in c(1:length(list_simul))){
    if(names(list_simul)[j]%in%l_name_2V){
      Simul_namev<-list_simul[[j]][,vect_times[Z]]
      Simul_ij[,Z]<-Simul_namev
      Obs_namev<-list_obs[[j]][,vect_times[Z]]
      Obs_ij[,Z]<-Obs_namev
      Z<-Z+1
    }
  }
  Obs_ij<-as.data.frame(Obs_ij)
  colnames(Obs_ij)<-sapply(c(1:length(l_name_2V)),function(x){
    return(paste0("V",x))
  })
  Simul_ij<-as.data.frame(Simul_ij)
  colnames(Simul_ij)<-colnames(Obs_ij)
  
  ### Empirical bivariate return levels
  Curves_found<-sapply(levels_used,FUN = function(x){
    texmex::JointExceedanceCurve(Sample = Obs_ij,
                                 ExceedanceProb = x)
  })
  Curves_found<-cbind.data.frame(Curves_found)
  All<-lapply(c(1:length(levels_used)),function(x){
    df_j<-data.frame(Curves_found[[x]])
    nj<-nrow(df_j)
    name_j<-rep(levels_used[x],nj)
    df_j$cl<-name_j
    return(df_j)
  })
  combined_df <- do.call(rbind, All)
  colnames(combined_df)<-c("t","s","level")
  combined_df$level<-as.character(combined_df$level)
  ### Bootstrap confidence regions
  Result_boot_confRegions<-replicate(500,
                                     Onesample_Bootstrap_Bivar_RL(Obs_t_s = Obs_ij,
                                                                  levels_used = levels_used))
  Df_transf<-apply(X=Result_boot_confRegions,
                   as.data.frame,MARGIN=2)
  m<-150
  vect_ind<-c(1:m)
  Theta_vector_target<-pi*(m+1-vect_ind)/(2*(m+1))
  Rult<-lapply(Theta_vector_target,Candidates_per_theta,
               list_boot_samples=Df_transf,
               hearts=hearts)
  Bounds_for_graph<-do.call(rbind.data.frame,Rult)
  lower <- Bounds_for_graph %>% filter(bounds == "Qinf")
  upper <- Bounds_for_graph %>% filter(bounds == "Qsup")
  
  # Build polygon data
  # x and y are changing so we must create a ggplot polygon object
  # arrange(level,t) for ordering
  # group_by level--> polygon per level
  # bind_rows to assemble rows per group
  poly_data <- lower %>%
    arrange(level, t) %>%
    group_by(level) %>%
    do({
      upper_part <- upper %>%
        filter(level == unique(.$level)) %>%
        arrange(desc(t))
      
      bind_rows(., upper_part)
    }) %>%
    ungroup()
  
  ### Curves for simulated
  Curves_found_SIM<-sapply(levels_used,FUN = function(x){
    texmex::JointExceedanceCurve(Sample = Simul_ij,
                                 ExceedanceProb = x)
  })
  Curves_found_SIM<-cbind.data.frame(Curves_found_SIM)
  All_sim<-lapply(c(1:length(levels_used)),function(x){
    df_j<-data.frame(Curves_found_SIM[[x]])
    nj<-nrow(df_j)
    ### Add levels used
    name_j<-rep(levels_used[x],nj)
    df_j$cl<-name_j
    return(df_j)
  })
  combined_df_sim<- do.call(rbind, All_sim)
  colnames(combined_df_sim)<-c("t","s","level")
  combined_df_sim$level<-as.character(combined_df_sim$level)
  ### Curves for data Texmex
  Texmex_curves<-list()
  d<-length(l_name_2V)
  Model_all <- mexAll(Obs_ij,mqu=0.7,dqu=rep(0.7,5))
  Simuls_from_tex_all <- mexMonteCarlo(nSample=5000,
                                       mexList=Model_all)
  Curves_tex_<-lapply(levels_used,FUN = function(x){
    JRC_j<-texmex::JointExceedanceCurve(
      Sample = Simuls_from_tex_all,
      ExceedanceProb = x,
      which=c("V1","V2"))
    
    return(JRC_j)
  })
  ### Create ggplot object
  GG_BivRL_obs_with_Tawn<-ggplot(Obs_ij)+
    geom_point(aes(x=V1,y=V2,col="data"),size=0.5)
  d_levels<-length(levels_used)
  linetypes <- scales::linetype_pal()(d_levels)
  for(j in c(1:d_levels)){
    GG_BivRL_obs_with_Tawn<-GG_BivRL_obs_with_Tawn+
      geom_jointExcCurve(x = Curves_tex_[[j]],aes(V1,V2,
                                                  col="CEVmodel"),
                         linetype=linetypes[j])
  }
  exp_type<-expression(x[M]^t~v)
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[1],
                                              x=l_name_2V[1],
                                              v=l_unit[[l_name_2V[1]]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[2],
                                              x=l_name_2V[2],
                                              v=l_unit[[l_name_2V[2]]])))
  GG_BivRL_obs<-GG_BivRL_obs_with_Tawn+
    geom_line(data=combined_df,aes(x=t,y=s,
                                   group=interaction(level),
                                   col="data",linetype=level))+
    geom_line(data=combined_df_sim,aes(x=t,y=s,
                                       group=interaction(level),
                                       col="simulations",
                                       linetype=level))+
    geom_polygon(data = poly_data,
                 aes(x = t,
                     y = s,
                     group = level),
                 alpha = 0.2,fill="lightgrey") +
    geom_line(data = Bounds_for_graph,
              aes(x = t,
                  y = s,
                  group = interaction(bounds, level),
                  linetype=level,col="confidence_band"))+
    labs(col="Legend",
         linetype="Level")+
    guides(fill="none")+
    xlab(xlab_title)+
    ylab(ylab_title)+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    scale_color_manual(values = cols_)
  return(GG_BivRL_obs)
}
Onesample_Bootstrap_Bivar_RL<-function(Obs_t_s,levels_used){
  Inds_t_s<-sample(x = c(1:nrow(Obs_t_s)),size = nrow(Obs_t_s),
                   replace = TRUE)
  New_data<-Obs_t_s[Inds_t_s,]
  Curves_found_1sample<-sapply(levels_used,FUN = function(x){
    texmex::JointExceedanceCurve(Sample =  New_data,
                                 ExceedanceProb = x)
  })
  Curves_found_1sample<-cbind.data.frame(
    Curves_found_1sample)
  All_1sample<-lapply(names(Curves_found_1sample),function(x){
    df_j<-data.frame(Curves_found_1sample[[x]])
    nj<-nrow(df_j)
    ### Add levels used
    name_j<-rep(x,nj)
    df_j$cl<-name_j
    return(df_j)
  })
  combined_df_1sample<- do.call(rbind, All_1sample)
  colnames(combined_df_1sample)<-c("t","s","level")
  combined_df_1sample$level<-levels_used[as.numeric(
    combined_df_1sample$level)]
  return(combined_df_1sample)
}


One_sample_EGPD<-function(n,theta_egpd,order_quantiles){
  kappa<-theta_egpd[["kappa"]]
  sig<-theta_egpd[["sigma"]]
  xi<-theta_egpd[["xi"]]
  simul_from_egpd<-mev::rextgp(n = n,type=1,kappa =kappa ,
                               sigma = sig,
                               xi =xi )
  model_boot_fit<-mev::fit.extgp(data=simul_from_egpd,
                                 model = 1,init = c(kappa,
                                                    sig,xi),
                                 method="mle",plots = FALSE)$fit$mle
  ### Quantiles found
  QEGPD<-sapply(order_quantiles,
                FUN=function(x){
                  return(mev::qextgp(p = x,
                                     kappa = model_boot_fit[["kappa"]],
                                     xi = model_boot_fit[["xi"]],
                                     sigma =model_boot_fit[["sigma"]]))
                })
  return(QEGPD)
}
Bootstrap_conf_band<-function(M,n,theta_egpd,order_quantiles,
                              alpha_){
  ### Estimaed EGPD curve
  Mean_egpd<-sapply(order_quantiles,
                    FUN=function(x){
                      return(mev::qextgp(p = x,
                                         kappa = theta_egpd[["kappa"]],
                                         xi = theta_egpd[["xi"]],
                                         sigma = theta_egpd[["sigma"]]))
                    })
  
  ### Conf_bands
  Result_boot<-t(replicate(n = M,
                           expr = One_sample_EGPD(n = n,
                                                  theta_egpd = theta_egpd,
                                                  order_quantiles = order_quantiles)))
  Niv_1<-alpha_/2
  Q1<-apply(
    Result_boot,MARGIN = 2,
    FUN = function(x){
      return(as.numeric(quantile(x,Niv_1)))
    })
  Niv_2<-1-(alpha_/2)
  Q2<-apply(Result_boot,MARGIN = 2,
            FUN = function(x){
              return(as.numeric(quantile(x,Niv_2)))
            })
  return(list("bound_inf"=Q1,
              "bound_sup"=Q2,
              "mean"=Mean_egpd))
}
Marg_2d_cond_ext<-function(list_simul,list_obs,vect_times,l_name_2V,
                           cols_ggplot,l_unit,risk_function,
                           ind_var_cond){
  
  Z<-1
  NSimul_for_compar<-nrow(list_simul[[1]])
  NObs_for_compar<-nrow(list_obs[[1]])
  Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
                   ncol =length(vect_times) )
  Obs_ij<-matrix(NA,nrow = NObs_for_compar,
                 ncol =length(vect_times) )
  for(j in c(1:length(list_simul))){
    if(names(list_simul)[j]%in%l_name_2V){
      Simul_namev<-list_simul[[j]][,vect_times[Z]]
      Simul_ij[,Z]<-Simul_namev
      Obs_namev<-list_obs[[j]][,vect_times[Z]]
      Obs_ij[,Z]<-Obs_namev
      Z<-Z+1
    }
  }
  Obs_ij<-as.data.frame(Obs_ij)
  colnames(Obs_ij)<-sapply(c(1:length(l_name_2V)),function(x){
    return(paste0("V_",x))
  })
  Simul_ij<-as.data.frame(Simul_ij)
  colnames(Simul_ij)<-colnames(Obs_ij)
  ### Select knowing that one coordinate is above the threshold
  Th_min<-min(apply(X = Obs_ij,MARGIN = 1,risk_function))
  TS_ind_ext<-which(Obs_ij[,ind_var_cond]>Th_min)
  TS_sim_ext<-which(Simul_ij[,ind_var_cond]>Th_min)
  ### Select sub dataset
  Sub_obs_ij<-Obs_ij[TS_ind_ext,]
  Sub_sim_ij<-Simul_ij[TS_sim_ext,]
  
  df_combo<-rbind.data.frame(Sub_obs_ij,Sub_sim_ij)
  df_combo$type<-c(rep("data",
                       nrow(Sub_obs_ij)),
                   rep("simulations",
                       nrow(Sub_sim_ij)))
  exp_type<-expression(x[M]^t~v)
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[1],
                                              x=l_name_2V[1],
                                              v=l_unit[[l_name_2V[1]]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t= vect_times[2],
                                              x=l_name_2V[2],
                                              v=l_unit[[l_name_2V[2]]])))
  GG0<-ggplot(data=Sub_obs_ij,aes(x=V_1,y=V_2,col="data"))+
    geom_point(size=2,pch=19)+
    geom_point(data=Sub_sim_ij,aes(x=V_1,y=V_2,col="simulations"),
               size=0.75,
               pch=17)+
    geom_xsidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    labs(col="Legend")+
    scale_color_manual(values=cols_ggplot)+
    scale_fill_manual(values= cols_ggplot)+
    guides(fill="none")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    xlab(xlab_title)+ylab(ylab_title)
  return(GG0)
  
}
# Extreme values estimators + data generated according to the model---------------------------------
#' Title
#'
#' @param gamma 
#' @param sigma 
#' @param kappa 
#' @param M 
#' @param order_quantiles 
#'
#' @return
#' @export
#'
#' @examples
RL_generations_EGPD<-function(gamma,sigma,kappa,M,order_quantiles,
                              list_params_EGPD){
  n.cyc<- list_params_EGPD[["n.cyc"]]
  mu.step<- list_params_EGPD[["mu.step"]]
  sigma.step<- list_params_EGPD[["sigma.step"]]
  nu.step<- list_params_EGPD[["nu.step"]]
  tau.step<- list_params_EGPD[["tau.step"]]
  Sample_EGPD<-mev::rextgp(n = M,xi= as.numeric(gamma),
                           kappa = kappa,sigma = sigma)
  Th<-quantile(Sample_EGPD,probs=0.50)
  #Th<-0
  ### PWM way
  SUB_EGPD<-Sample_EGPD[which(Sample_EGPD>Th)]-Th
  mu0<-mean(SUB_EGPD)
  next_<-length(SUB_EGPD)+1
  Fbar_side<-(next_-VineCopula::pobs(SUB_EGPD)*next_)/next_
  mu1<-mean(SUB_EGPD*Fbar_side)
  Shape<-(mu0-4*mu1)/(mu0-2*mu1)
  print(Shape)
  Scale<-mu0*(1-Shape)
  # FIT_pos<-mev::gp.fit(xdat = Sample_EGPD,threshold =Th)$est
  # Shape<-FIT_pos[2]
  # #Use gpd property to get the scale at 0.
  # Scale<-FIT_pos[1]-Th*Shape
  # #Scale<-FIT_pos[1]
  # INIT<-c(Shape,Scale)
  #Nu.start initialisation using moments.
  Th_beg<-quantile(x = Sample_EGPD,0.10)
  Lower_tail<-Sample_EGPD[which(Sample_EGPD<Th_beg)]
  Moment_1<-mean(Lower_tail)
  nu.start<-as.numeric((1-(Moment_1/Th_beg))^(-1)-1)
  db<-as.data.frame(Sample_EGPD)
  colnames(db)<-c("x")
  con <- gamlss::gamlss.control(n.cyc =  n.cyc,
                                mu.step = mu.step, sigma.step = sigma.step, 
                                nu.step = nu.step,tau.step = tau.step,autostep=TRUE,
                                trace = TRUE)
  con.i<-gamlss::glim.control(glm.trace = TRUE)
  EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")
  Fitting_sample_EGPD <- gamlss::gamlss(x~1, 
                                        data=db,
                                        family = EGPD1Family(mu.link = "identity"),
                                        control = con,mu.start=gamma,
                                        sigma.start=sigma,
                                        nu.start=kappa,
                                        i.control=con.i,
                                        method=CG())
  muFit <- fitted(Fitting_sample_EGPD,"mu")[1]
  sigmaFit <- predict(Fitting_sample_EGPD,what="sigma", 
                      type="response")[[1]]
  nuFit <- predict(Fitting_sample_EGPD,what="nu", 
                   type="response")[[1]]
  # Ratio_found_EGPD<-Ratio_log_model1(gamma = muFit,sigma_EGPD = sigmaFit,order_quantiles = order_quantiles,
  #                  data = Sample_EGPD)
  RL_EGPD<-sapply(order_quantiles,function(x){
    return(as.numeric(mev::qextgp(p =x,
                                  kappa = nuFit,sigma = sigmaFit,xi =muFit)))
  })
  return(RL_EGPD)
}
One_sim<-function(ALL_Quantiles,Kappa_t,gam_t,sig_t,list_params_EGPD,
                  Msim,i){
  return(RL_generations_EGPD(kappa=Kappa_t,
                             gamma=gam_t,sigma=sig_t,
                             list_params_EGPD= list_params_EGPD,
                             M=Msim,order_quantiles = ALL_Quantiles))
}
# Extreme values estimators + data generated according to the model---------------------------------
#' Title
#'
#' @param gamma 
#' @param sigma 
#' @param kappa 
#' @param M 
#' @param order_quantiles 
#'
#' @return
#' @export
#'
#' @examples
RLevel_generations_EGPD<-function(gamma,sigma,kappa,M,order_quantiles){
  
  Sample_EGPD<-rEGPDModel1(n = M,mu = gamma,nu = kappa,sigma = sigma)
  Theta<-mev::gp.fit(xdat = Sample_EGPD,threshold = 0)$est
  INIT<-c(Theta[2],Theta[1])
  db<-as.data.frame(Sample_EGPD)
  colnames(db)<-c("x")
  con <- gamlss.control(n.cyc = 100,
                        mu.step = 0.1, 
                        sigma.step = 0.1, nu.step = 0.1,tau.step = 0.1,autostep=TRUE,
                        trace = FALSE)
  con.i=glim.control(glm.trace = FALSE)
  Fitting_sample_EGPD <- gamlss(x~1, 
                                data=db,
                                family = EGPD1Family(mu.link = "identity"),
                                control = con,mu.start=INIT[1],sigma.start=INIT[2],
                                nu.start=0.5,
                                i.control=con.i,
                                method=CG())
  muFit <- fitted(Fitting_sample_EGPD,"mu")[1]
  sigmaFit <- predict(Fitting_sample_EGPD,what="sigma", 
                      type="response")[[1]]
  nuFit <- predict(Fitting_sample_EGPD,what="nu", 
                   type="response")[[1]]
  Rl_estim<-RL_EGPD_model1(gamma = muFit,sigma_EGPD = sigmaFit,
                           order_quantiles = order_quantiles,
                           Kappa = nuFit)
  return(Rl_estim)
}
### 
#' Gaussian_mixture vs dnormal distrib
#'
#' @param vect_Id: matrix. Simulations from a N(0,1)
#' @param alpha: float. Confidence level
#'
#' @return Indicator of null hypothesis rejection
#' @export
#'
#' @examples
TwoGM_versus_normal<-function(vect_Id,alpha){
  phi1_mid<-Phi1_normal_stat(vect_Id = vect_Id,
                             alpha = alpha/2)
  phi2_mid<-Phi2_normal_stat(vect_Id = vect_Id,
                             alpha = alpha/2)
  l_2_based<-max(phi1_mid,phi2_mid)
  phi4_mid<-Phi4_normal_stat(vect_Id = vect_Id,
                             alpha = alpha/2)
  phi5_mid<-Phi5_normal_stat(vect_Id = vect_Id,
                             alpha = alpha/2)
  l_inf_based<-max(phi4_mid,phi5_mid)
  return(list("l_inf"=l_inf_based,
              "l_2"=l_2_based))
  
}
Phi1_normal_stat<-function(vect_Id,alpha){
  n<-nrow(vect_Id)
  d<-ncol(vect_Id)
  mu_data<-colMeans(vect_Id)
  S<-t(mu_data)%*%mu_data*n
  Thresh_null<-qchisq(p = 1-(alpha),df = d)
  return(as.numeric(S>Thresh_null))
}
Phi2_normal_stat<-function(vect_Id,alpha){
  n<-nrow(vect_Id)
  d<-ncol(vect_Id)
  S<-rowSums(vect_Id**2)
  Thresh_null<-qchisq(p=1-(alpha/n),df=d)
  return(max(as.numeric(S>Thresh_null)))
}
### l(inf) based
Phi4_normal_stat<-function(vect_Id,alpha){
  n<-nrow(vect_Id)
  d<-ncol(vect_Id)
  S<-n*colMeans(vect_Id)**2
  Thresh_null<-qchisq(p=1-(alpha/d),df=1)
  return(max(as.numeric(S>Thresh_null)))
}
Phi5_normal_stat<-function(vect_Id,alpha){
  n<-nrow(vect_Id)
  d<-ncol(vect_Id)
  Ind_Chi_Hnull<-vect_Id**2
  Thresh_null<-qchisq(p=1-(alpha/(d*n)),df=1)
  R<-Ind_Chi_Hnull>Thresh_null
  Max_per_j<-apply(R,MARGIN = 2,FUN = max)
  return(max(Max_per_j))
}
### Alternative hypothesis
Simul_Contamination_Model<-function(prob_classes,nu){
  d<-length(nu)
  cl_used<-sample(c(1:length(prob_classes)),prob = prob_classes,
                  size = 1)
  if(cl_used==1){
    u<-mvtnorm::rmvnorm(n = 1,mean = rep(0,d))
  }else{
    u<-mvtnorm::rmvnorm(n = 1,mean = nu)
  }
  return(u)
}
Candidates_per_theta<-function(list_boot_samples,Theta,alpha_prop,
                               hearts){
  Rult_for_theta<-parLapply(cl=hearts,list_boot_samples,
                            fun = Candidates_per_theta_1sample,
                            Theta=Theta)
  Rult_for_theta<-do.call(rbind.data.frame,Rult_for_theta)
  Rult_for_theta$level<-as.character(Rult_for_theta$level)
  Rult_for_theta$R<-rowSums(cbind(Rult_for_theta$t,
                                  Rult_for_theta$s)**2)**(1/2)
  Rays_found<-Rult_for_theta %>% 
    group_by(level) %>% 
    summarise(Qinf=quantile(R,0.025),
              Qsup=quantile(R,0.975))
  X<-Rays_found[,c(2:3)]*cos(Theta)
  Y<-Rays_found[,c(2:3)]*sin(Theta)
  Levels<-Rays_found[,1]
  
  Combinaison_coords<-cbind(melt(X,
                                 id.vars = NULL),
                            melt(Y,
                                 id.vars = NULL)[,2])
  colnames(Combinaison_coords)<-c("bounds","t","s")
  ncomb<-nrow(Combinaison_coords)
  n_for_rep<-ncomb/nrow(X)
  Combinaison_coords$level<-rep(unlist(Levels),n_for_rep)
  Combinaison_coords$theta<-rep(Theta,ncomb)
  return(Combinaison_coords)
}

Candidates_per_theta_1sample<-function(Boot_Sample,Theta){
  Boot_Sample$angle<-atan(Boot_Sample[,2]/Boot_Sample[,1])
  Inds_chosen<-Boot_Sample %>% group_by(level) %>%
    summarise(t=t[which.min(abs(Theta-angle))],
              s=s[which.min(abs(Theta-angle))])
  return(sapply(X =Inds_chosen,FUN = unlist))
}
Params_HTawn_one_dqu<-function(mqu,dqu,vect_l,ind_ref){
  model_texmex<-texmex::mex(vect_l,mqu = mqu,
                            dqu = dqu,which = ind_ref)
  Dep<-model_texmex$dependence
  Theta<-Dep$coefficients
  Z<-Dep$Z
  n<-nrow(Z)
  d<-ncol(Theta)
  Ref_var<- seq(dqu, 1 - 1/n, length = n)
  Result_per_column<-NA
  Result_per_column2<-NA
  if(d>1){
    ### Verify if test doable
    Sum_ind<-sum(as.numeric(!is.na(Theta[1,])))
    n<-nrow(Z)
    if(Sum_ind==d){
      Resid<-Z
      Absolute_values<-abs(Resid-colMeans(Resid))
      Result_per_column<-apply(Resid,MARGIN = 2,
                               FUN = function(x){
                                 result<-round(cor.test(x,Ref_var,method="kendall")$p.value,
                                               3)
                                 return(result)
                               })
      Result_per_column2<-apply(X = Absolute_values,MARGIN = 2,
                                FUN = function(x){
                                  result<-round(cor.test(x,Ref_var,method="kendall")$p.value,
                                                3)
                                  return(result)
                                })
    }
  }else{
    ### Verify if test doable
    Answer<-!is.na(as.numeric(Theta[1,]))
    if(Answer==TRUE){
      Resid<-as.numeric(Z)
      Absolute_values<-abs(Resid-mean(Resid))
      Result_per_column<-round(cor.test(Resid,
                                        Ref_var, method="kendall")$p.value,3)
      Result_per_column2<-round(cor.test(Absolute_values,
                                         Ref_var,method="kendall")$p.value,3)
    }
  }
  return(list("Test_1"=Result_per_column,
              "Test_2"=Result_per_column2,
              "Theta"=Theta))
}
Analysis_diag_HTawn_Conv_Hull<-function(Theta_opt){
  d<-length(Theta_opt)
  list_GG_for_hull<-list()
  for(Z in c(1:d)){
    HTawn_z<-Theta_opt[[Z]]
    Modelboot_z<-texmex::bootmex(x = HTawn_z)
    ### Boot run
    j<-1
    Results<-t(sapply(Modelboot_z$boot,function(x){
      y<-x$dependence
      return(c(y[1,j],y[2,j]))
    }))
    list_GG_for_hull[[Z]]<-as.data.frame(Results)
  }
  ### Convex graph
  Combinaison_for_hull<-as.data.frame(
    matrix(t(data.frame(list_GG_for_hull)),
           ncol = 2))
  B_rep<-nrow(Combinaison_for_hull)/length(l_name)
  
  colnames(Combinaison_for_hull)<-c("a","b")
  Combinaison_for_hull$variable_cond<-as.character(
    c(sapply(X = c(1:length(l_name)),
             FUN = function(x){return(rep(x,B_rep))}))
  )
  Chull<-Combinaison_for_hull %>% 
    group_by(variable_cond) %>% 
    slice(chull(a, b))
  
  GG_hull<-ggplot(Combinaison_for_hull, 
                  aes(a, b)) +
    geom_point() +
    geom_polygon(data =Chull, alpha = 0.3,aes(fill=variable_cond)) +
    theme_minimal()+labs(fill="Legend")+
    theme(axis.title=element_text(size=25),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          axis.text = element_text(size=12),
          strip.text = element_text(size = 12))
  return(GG_hull)
}
Analysis_diag_HTawn_evol_DQU<-function(vect_dqu,
                                       Vect_obs,MQU,
                                       chosen_dqu){
  d<-ncol(Vect_obs)
  GG_theta<-ggplot()
  Rult_simplified<-c()
  Indep_whole<-c()
  for(cond_var in c(1:d)){
    Results_Indep_Theta<-sapply(vect_dqu,
                                Params_HTawn_one_dqu,
                                mqu=MQU,vect_l=Vect_obs,
                                ind_ref=cond_var)
    Inds_for_analyse<-apply(X = Results_Indep_Theta,
                            FUN = function(x){
                              Ind_nafalse<-sum(rowSums(!is.na(x$Theta)))
                              dim<-ncol(x$Theta)*nrow(x$Theta)
                              return(Ind_nafalse==dim)
                            },MARGIN = 2)
    ThetaI_simplified<-Results_Indep_Theta[,Inds_for_analyse]
    Sub_DQU<-vect_dqu[Inds_for_analyse]
    Result_Independence_test<-t(ThetaI_simplified[c(1:2),])
    Result_several_theta<-do.call(rbind,ThetaI_simplified[3,])
    ### Drop NA results
    
    Rult<-melt(Result_several_theta)
    Sub_vars<-c("a","b")
    Nb_rep_DQU<-length(Sub_vars)
    Inds_chosen<-which(Rult$Var1%in%Sub_vars)
    Rult2<-Rult[Inds_chosen,]
    N_rep<-nrow(Rult2)/(Nb_rep_DQU*length(Sub_DQU))
    colnames(Rult2)<-c("param","variable","value")
    chosen_dquj<-chosen_dqu[cond_var]
    Rult2$dqu<-rep(chosen_dquj,
                   nrow(Rult2))
    Variable_used<-Rult2$variable[1]
    Rult2$quantile<-rep(c(sapply(Sub_DQU,
                                 FUN = function(x){rep(x,Nb_rep_DQU)})),
                        N_rep)
    Rult_simplified<-rbind(Rult_simplified,
                           Rult2)
    
    ### Independence test
    ######################
    Melting_indep<-melt(apply(Result_Independence_test,
                              MARGIN = 2,FUN = unlist))
    colnames(Melting_indep)<-c("variable","test_used",
                               "pval")
    Melting_indep$test_used<-ifelse(Melting_indep$test_used=="Test_1",
                                    yes = "Z vs Y",
                                    no = "|Z-mean(Z)| vs Y")
    Melting_indep$quantile<-rep(c(sapply(Sub_DQU,
                                         FUN = function(x){rep(x,N_rep)})),
                                Nb_rep_DQU)
    Melting_indep$variable_cond<-rep(Variable_used,
                                     nrow(Melting_indep))
    Melting_indep$dqu<-rep(chosen_dquj,
                           nrow(Melting_indep))
    Indep_whole<-rbind(Indep_whole,
                       Melting_indep)
  }
  GG_theta<-ggplot(Rult_simplified
                   ,aes(x=quantile,y=value,
                        group=interaction(variable,param),
                        col=param))+
    geom_line()+
    ylab("Estimator")+
    xlab("Dqu")+
    facet_wrap(~variable,
               scales = "free_y")+
    geom_point()+
    labs(col="Legend")+
    geom_vline(aes(xintercept=dqu))
  GG_Indep<-ggplot(Indep_whole,
                   aes(x=quantile,y=pval,
                       group=interaction(variable_cond,
                                         test_used),
                       col=variable_cond))+
    geom_line()+
    facet_wrap(~test_used)+
    geom_point()+
    ylab("p value")+
    xlab("Dqu")+
    labs(col="Legend")+
    geom_hline(yintercept = 0.05,col="black")+
    geom_vline(aes(xintercept = dqu,
                   col=variable_cond))
  if(d==1){
    GG_Indep<-GG_Indep+
      guides(col="none")
  }
  return(list("theta"=GG_theta,
              "indep"=GG_Indep))
}
