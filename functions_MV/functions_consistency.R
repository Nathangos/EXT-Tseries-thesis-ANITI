
#' Evolution with the time lag 
#' of the correlation coefficient for each pair of variables.
#'
#' @param list_simul list[str: dataframe]. Simulated multivariate
#' time series where each value corresponds to the 
#' measures of one variable. 
#' @param list_reality list[str: dataframe]. Observed multivariate
#' time series where each value corresponds to the 
#' measures of one variable. 
#' @param q_per_time vector[float]. Proportion of 
#' extreme values.
#' @param B int. Number of bootstrap samples to provide
#' a confidence band
#' @param list_names_VAR vector[string]. Variable names
#'
#' @return Dataframe. 
#' @export
#'
#' @examples
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
   
  },MARGIN = 1)
  return(All_C_extremo)
}
#' Evolution with the time lag 
#' of the correlation coefficient for a given pair (name_cond, name_other)
#'
#' @param list_simul list[str: dataframe]. Simulated multivariate
#' time series where each value corresponds to the 
#' measures of one variable. 
#' @param list_reality list[str: dataframe]. Observed multivariate
#' time series where each value corresponds to the 
#' measures of one variable. 
#' @param name_cond string. First variable. 
#' @param name_other string. Second variable. 
#' @param q_per_time vector[float]. Proportion of 
#' extreme values.
#' @param B int. Number of bootstrap samples to provide
#' a confidence band
#'
#' @return Dataframe. 
#' @export
#'
#' @examples
Cross_Extremo_MV_ij<-function(list_simul,list_reality,name_cond,name_other,
                           q_per_time,B){
  
  ### Obs
  Variable_reality_cond<-list_reality[[name_cond]]
  Variable_reality_other<-list_reality[[name_other]]
  
  ### Simul
  Variable_simul_cond<-list_simul[[name_cond]]
  Variable_simul_other<-list_simul[[name_other]]
  
  ### Modify according to the available times.  
  Intersectg_times<-intersect(colnames(Variable_reality_cond),
                              colnames(Variable_reality_other))
  L<-length(Intersectg_times)
  Variable_reality_cond<-Variable_reality_cond[,Intersectg_times]
  colnames(Variable_reality_cond)<-c(1:L)
  Variable_reality_other<-Variable_reality_other[,Intersectg_times]
  colnames(Variable_reality_other)<-colnames(Variable_reality_cond)
  ##
  Variable_simul_cond<-Variable_simul_cond[,Intersectg_times]
  colnames(Variable_simul_cond)<-colnames(Variable_reality_cond)
  Variable_simul_other<-Variable_simul_other[,Intersectg_times]
  colnames(Variable_simul_other)<-colnames(Variable_reality_cond)
  z<-1
  vect_distances<-c()
  Matrix_pairs<-list()
  
  for(j in 1:L){
    valeurs_t_plus_h<-j:L
    for (i in valeurs_t_plus_h){
      Matrix_pairs[[z]]<-c(i,j)
      vect_distances<-c(vect_distances,abs(i-j))
      z<-z+1
    }
  }
  Tau_cond<-sapply(c(1:L),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_reality_cond)
  Tau_other<-sapply(c(1:L),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_reality_other)
  list_Tau<-list("other"=Tau_other,"cond"=Tau_cond)
  ##
  CROSS_extremo_reality<-cross_extremogram(Matrix_pairs = Matrix_pairs,
                         var_cond = Variable_reality_cond,
                         var_other = Variable_reality_other,
                         l_Tau = list_Tau)
  ### Attention ! We use then the quantile from simulations
  Tau_cond_sim<-sapply(c(1:L),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_per_time,
    x=Variable_simul_cond)
  Tau_other_sim<-sapply(c(1:L),function(x,q,t){
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

#' Value of the cross-extremogram for a given 
#' pair and resampled observations. 
#'
#' @param vect_distances vector[float]. Distance 
#' between time steps.
#' @param Matrix_pairs List. Each possible
#' time pairs.
#' @param var_cond dataframe. Univariate time series
#' for the first variable.
#' @param var_other dataframe. Univariate time series
#' for the other variable.
#' @param l_Tau list. It contains for each time series
#' the vector of threshold values. 
#'
#' @return Dataframe with the evolution of the cross-extremogram
#' with the time lag. 
#' @export
#'
#' @examples
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
#' Obtain from a n-sample its (1-Alpha_q)% bounds.  
#'
#' @param x vector[float]. Observed values. 
#' @param Alpha_q float. Confidence level.
#'
#' @return vector[float]. Bounds of the confidence interval
#' @export
#'
#' @examples
F_bounds_per_column<-function(x,Alpha_q){
  
  return(c(as.numeric(
    quantile(x,
             c(Alpha_q/2,1-(Alpha_q/2)
             )))))
}
#' Evolution of the dependence structure between two time steps 
#' with the time lag 
#' 
#' @param Matrix_pairs List. Each possible
#' time pairs.
#' @param var_cond Dataframe. Univariate time series 
#' for one variable.
#' @param var_other Dataframe. Univariate time series 
#' for the other variable.
#' @param l_Tau list. It contains for each time series
#' the vector of threshold values. 
#'
#' @return
#' @export
#'
#' @examples
cross_extremogram<-function(Matrix_pairs,var_cond,var_other,l_Tau){
  
  PIST_empirique<-sapply(Matrix_pairs,FUN = cross_pi_s_t,
                         var_cond=var_cond,
                         var_other=var_other,l_Tau=l_Tau)
  return(PIST_empirique)
}
#' Extremal correlation coefficient between 
#' two time steps, for two specific components. 
#'
#' @param couple_s_t vector. Pair of time steps. 
#' @param var_cond dataframe. Univariate time series
#' for the first variable.
#' @param var_other dataframe. Univariate time series
#' for the other variable.
#' @param l_Tau list. It contains for each time series
#' the vector of threshold values. 
#'
#' @return float. Value of the empirical estimator.
#' @export
#'
#' @examples
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
### Work on the distribution of angular time series -----------
################

#' F_automatic_ggplot_marg_Omega
#'
#' @param Obj_data: df. Concatenated time series. 
#' @param Obj_sim: df. Simulated concatenated angles.
#' @param l_name_time: list[str:int]. List containing
#' for each selected variable the analysed time. 
#' @param cols_ggplot: vector[str]. Vector of desired cols in 
#' the ggplot plot. 
#' @param list_CumSum_Time: list[float]. List associating
#' each forcing condition to the cumulated number of times observed
#'
#' @return GGplot comparing the bivariate distribution of simulated angle 
#' and observed from extreme time series. 
#' @export
#'
#' @examples
F_automatic_ggplot_marg_Omega<-function(Obj_data,
                        Obj_sim,l_name_time,cols_ggplot,
                        list_CumSum_Time,l_unit){
  
  name_first<-names(l_name_time)[1]
  name_second<-names(l_name_time)[2]
  tprime<-l_name_time[[1]]
  sprime<-l_name_time[[2]]
  ### Select the corresponding time in the 
  ### concatenated angles using the number of 
  ### obs times for each forcing condition

  DF1_analysed_sim<-Obj_sim[[name_first]]
  DF2_analysed_sim<-Obj_sim[[name_second]]
  d1<-ncol(DF1_analysed_sim)
  d2<-ncol(DF2_analysed_sim)
  if((!tprime%in%c(1:d1))|(!sprime%in%c(1:d2))){
    return(NA)
  }
  U_sim<-DF1_analysed_sim[,tprime]
  V_sim<-DF2_analysed_sim[,sprime]
  
  DF1_analysed_data<-Obj_data[[name_first]]
  DF2_analysed_data<-Obj_data[[name_second]]
  U<-DF1_analysed_data[,tprime]
  V<-DF2_analysed_data[,sprime]
  DF_1<-cbind.data.frame(U,V)
  DF_sim<-cbind.data.frame(U_sim,V_sim)
  colnames(DF_1)<-sapply(c(1:2),function(x){return(paste0("V",x))})
  colnames(DF_sim)<-colnames(DF_1)
  df_combo<-rbind.data.frame(DF_1,DF_sim)
  df_combo$Legend<-c(rep("data",nrow(DF_1)),
                     rep("simulations",nrow(DF_sim)))
  
  
  exp_type<-expression(x[M]^t~v)
  name_first_correct<-Fct_correct_name(x = name_first,
                             target = "Surcote",
                             replacement = "Surge")
  name_second_correct<-Fct_correct_name(x = name_second,
                              target = "Surcote",
                              replacement = "Surge")
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                         list(t=tprime,
                              x=name_first,
                              v=l_unit[[name_first]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                         list(t=sprime,
                              x=name_second_correct,
                              v=l_unit[[name_second]])))
  Name_file<-paste0("/",CPLMT_graph_theta,"_",METHOD_ANGLE,
                    "_Distrib_Omega_",name_first,"_t=",
                    tprime,"_",name_second,"_s=",sprime,"_d=",
                    length(l_unit),".png")
  GG1<-ggplot2::ggplot(data=df_combo,aes(x=V1,y=V2,colour=Legend,shape=Legend),
  )+
    geom_point(aes(size=Legend))+
    geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    
    scale_color_manual(values=cols_ggplot)+
    scale_fill_manual(values=cols_ggplot)+
    scale_shape_manual(values = c("simulations"=17,"data"=19))+
    scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
    guides(fill="none")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    xlab(xlab_title)+ylab(ylab_title)
  ggsave(filename = paste0(Theta_path,
                           Name_file),
         plot = GG1,
         width=8,height=6)
}
#' Comparing graphically observed and simulated time
#' series using several tools. 
#'
#' @param l_name_time list. Variable and 
#' time step analysed.
#' @param l_unit vector[string]. Unit of each variable.
#' @param risk_function function. Risk function used 
#' in one of the tools. 
#' @param ind_var_cond int. Conditioning variable
#' used in one of the ggplot objects. 
#' @param levels_used Float. Level of the return curve.
#' @param cols_ vector [float]. 
#' @param L_simul list. Simulated extreme multivariate time series. 
#' @param L_obs list. Observed extreme multivariate time series. 
#' @param L_whole_obs list. Observed multivariate time series (standard and extreme). 
#' @param Path_TREND String. Directory path for exporting ggplot objects.
#' @param Path_EXT String. Directory path for exporting ggplot objects.
#' @param hearts String. Other directory path for exporting ggplot objects.
#' @param basis_functionGAM 
#' @param n_sim int. Number of bootstrapped sampled 
#' 
#' @param list_info_export 
#' @param theme_common Theme object from ggplot2 library.
#'
#' @return NA
#' @export
#'
#' @examples
F_automatic_ggplot_marg<-function(l_name_time, l_unit,
            risk_function,ind_var_cond,levels_used,cols_,
            L_simul,L_obs,L_whole_obs,
            Path_TREND,Path_EXT,hearts,basis_functionGAM,
            n_sim,list_info_export,theme_common){
  
  Prefix_export<-list_info_export[["prefix"]]
  CPLMT_export<-list_info_export[["complement"]]
  opt_used<-list_info_export[["model"]]
  Z<-1
  NSimul_for_compar<-nrow(L_simul[[1]])
  NObs_for_compar<-nrow(L_obs[[1]])
  Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
                   ncol = length(l_name_time) )
  Obs_ij<-matrix(NA,nrow = NObs_for_compar,
                 ncol = length(l_name_time) )
  NAMES<-names(l_name_time)
  D<-length(L_simul)
  D_name_time<-length(l_name_time)
  for(j in c(1:D_name_time)){
    
    Key_t<-NAMES[j]
    if(Key_t%in%names(L_simul)){
      value_t<-l_name_time[[j]]
      Simul_namev<-L_simul[[Key_t]][,as.character(value_t)]
      Simul_ij[,j]<-Simul_namev
      Obs_namev<-L_obs[[Key_t]][,as.character(value_t)]
      Obs_ij[,j]<-Obs_namev
    }
  }
  Obs_ij<-as.data.frame(Obs_ij)
  colnames(Obs_ij)<-sapply(c(1:length(l_name_time)),function(x){
    return(paste0("V_",x))
  })
  Simul_ij<-as.data.frame(Simul_ij)
  colnames(Simul_ij)<-colnames(Obs_ij)
  Names_vars<-names(l_name_time)
  T1<-as.character(l_name_time[[1]])
  T2<-as.character(l_name_time[[2]])
  N1<-Names_vars[1]
  N2<-Names_vars[2]
  obj<-Marg_2d_simul_vs_obs(l_name_time = l_name_time,
                  Obs_ij = Obs_ij,Simul_ij = Simul_ij,
                  cols_ggplot = cols_,l_unit = l_unit)+
    theme_common
  F_name1<-paste0(Path_TREND
                ,"/",Prefix_export,"_distrib_values_",N1,"t=",
                T1,"X",N2,"s=",
                T2,CPLMT_export,"_d=",D,
                opt_used,".png")
  ggsave(filename = F_name1,
         width=8,height = 6,
         plot = obj)
  
  ## Bivariate return level
  ##### on the original data observations.
  N_ext<-nrow(L_obs[[1]])
  N_full<-nrow(L_whole_obs[[1]])
  zeta_forsim<-N_ext/N_full
  Obs_with_cond<-cbind.data.frame(L_obs[[N1]][,T1],
                                  L_obs[[N2]][,T2])
  colnames(Obs_with_cond)<-sapply(c(1:length(l_name_time)),function(x){
    return(paste0("V",x))
  })
  Obs_without_cond<-cbind.data.frame(L_whole_obs[[N1]][,T1],
                                     L_whole_obs[[N2]][,T2])
  Obs_without_cond<-as.data.frame(Obs_without_cond)
  colnames(Obs_without_cond)<-sapply(c(1:length(l_name_time)),function(x){
    return(paste0("V",x))
  })
  BivRL<-Bivariate_RLevel_simul_vs_obs(Obs_ij = Obs_without_cond,
           Simul_ij = Simul_ij, l_name_time = l_name_time,
           cols_ggplot = cols_,zeta_sim = zeta_forsim,
           l_unit = l_unit,levels_used = levels_used,
           hearts = hearts,Obs_ext_ij = Obs_with_cond,
           Nb_boot = 500)+
    theme_common
  F_name2<-paste0(Path_EXT
                  ,"/",Prefix_export,"BivReturnLevel_",N1,"t=",
                  T1,"X",N2,"t=",
                  T2,CPLMT_export,"_d=",D,
                  opt_used,".png")
  ggsave(filename = F_name2,
         width=10,height = 6,
         plot = BivRL)
  
  ### Threshold for bivariate quantile regressions
  # Seq_th<-seq.int(from = 0.5,to = 0.95,length.out =2)
  # list_ggplot_qreg<-Bivariate_Qreg_simul_vs_obs(
  #         Obs_ij = Obs_without_cond,
  #         Simul_ij = Simul_ij,l_name_time=l_name_time,
  #         cols_ggplot = cols_,l_unit = l_unit,
  #          zeta_sim = zeta_forsim,basis_functionGAM = basis_functionGAM,
  #          hearts = hearts,n_sim=n_sim,Seq_th=Seq_th)
  # L_gg_plot<-length(list_ggplot_qreg)
  # for(j in c(1:L_gg_plot)){
  #   plot_elt<-list_ggplot_qreg[[j]]+
  #     theme_common
  #   F_namej<-paste0(Path_EXT,"/",Prefix_export,
  #                   "qreg_",N1,"t=",
  #                   T1,"X",N2,"t=",
  #                   T2,CPLMT_export,"_d=",D,
  #                   opt_used,"_",j,".png")
  #   ggsave(filename = F_namej,
  #          width=8,height = 6,
  #          plot = plot_elt)
  # }
  # 
  ## Cond_ext
  obj_cond<-Marg_2d_cond_ext(Simul_ij = Simul_ij,
                             Obs_ij = Obs_ij,
                             l_name_time = l_name_time ,
                             cols_ggplot=cols_,
                             l_unit = l_unit,
                             risk_function =risk_function,
                             ind_var_cond = ind_var_cond)+
    theme_common
  F_name_END<-paste0(Path_TREND,"/",Prefix_export,
                     "cond_distrib_values_",N1,"t=",
                     T1,"X",N2,"t=",
                     T2,CPLMT_export,"_d=",D,
                     opt_used,".png")
  ggsave(filename = F_name_END,
         width=8,height = 6,
         plot = obj_cond)
}

#' Compare graphically the simulated and observed
#' extreme time series for a given pair of time steps.
#'
#' @param l_name_time list. Variable and 
#' time step analysed.
#' @param Obs_ij dataframe. Observed values at the given pair.
#' @param Simul_ij dataframe. Simulated values at the given pair. 
#' @param cols_ggplot vector[string]. Graphical options
#' to distinguish observations and simulations.
#' @param l_unit vector[string]. Unit of each variable.
#'
#' @return GGplot.
#' @export
#'
#' @examples
Marg_2d_simul_vs_obs<-function(l_name_time,Obs_ij,Simul_ij,
                               cols_ggplot,l_unit){
  
  
  NAMES<-names(l_name_time)
  df_combo<-rbind.data.frame(Obs_ij,Simul_ij)
  df_combo$type<-c(rep("data",nrow(Obs_ij)),
                   rep("simulations",nrow(Simul_ij)))
  exp_type<-expression(x[M]^t~v)
  N1<-NAMES[1]
  N2<-NAMES[2]
  T1<-l_name_time[[N1]]
  T2<-l_name_time[[N2]]
  Vect_corrected<-sapply(c(N1,N2),
                  FUN = function(x){
            Fct_correct_name(x = x,
                        target = "Surcote",
                        replacement = "Surge")
                  })
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                    list(t=T1,
                         x=Vect_corrected[1],
                         v=l_unit[[N1]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                       list(t=T2,
                            x=Vect_corrected[2],
                            v=l_unit[[N2]])))
  GG0<-ggplot2::ggplot(data=Obs_ij,aes(x=V_1,y=V_2,col="data"))+
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
Bivariate_Qreg_simul_vs_obs<-function(l_name_time,Obs_ij,
                                      Simul_ij,cols_ggplot,l_unit,
                                      Seq_th,zeta_sim,
                                      basis_functionGAM,
                                      hearts,n_sim){

  Obs_ij<-as.data.frame(Obs_ij)
  colnames(Obs_ij)<-c("X_regressor","Y_target")
  
  Simul_ij<-as.data.frame(Simul_ij)
  colnames(Simul_ij)<-c("X_regressor","Y_target")
  
  ### Sample from observations regressors
  Ind_data_new<-sample(c(1:nrow(Obs_ij)),
           size=n_sim,replace = TRUE)
  Data_new<-Obs_ij[Ind_data_new,]
  
  ### Observations
  Df_ald_qregs<-lapply(Seq_th,ALD_reg_XY,
              data_XY = Obs_ij,
             basis_functionGAM = basis_functionGAM,
             Sim_X = Data_new)
  Combinations_orig<-do.call(rbind.data.frame,
                     args = Df_ald_qregs)
  ### Simulations
  ### Sample from simulations regressors
  Ind_simul_new<-sample(c(1:nrow(Simul_ij)),
                        size=n_sim,replace = TRUE)
  Sim_new<-Simul_ij[Ind_simul_new,]
  Df_ald_qregs_sims<-lapply(Seq_th,ALD_reg_XY,
              data_XY = Simul_ij,
              basis_functionGAM = basis_functionGAM,
              Sim_X = Sim_new)
  Combinations_fr_simul<-do.call(rbind.data.frame,
                              args = Df_ald_qregs_sims)
  Combinations<-rbind(Combinations_orig,
                      Combinations_fr_simul)
  Combinations$group<-c(rep("data",
                            nrow(Combinations_orig)),
        rep("simulations",nrow(Combinations_fr_simul)))
 
  All_linetypes<-c("solid", "dashed", "dotted",
              "dotdash", "longdash", "twodash")
  exp_type<-expression(x[M]^t~v)
  NAMES<-names(l_name_time)
  N1<-NAMES[1]
  N2<-NAMES[2]
  T1<-l_name_time[[N1]]
  T2<-l_name_time[[N2]]
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t=T1,
                                              x=N1,v=l_unit[[N1]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t=T2,
                                              x=N2,v=l_unit[[N2]])))
  list_gg<-list()
  for(j in c(1:length(unique(Combinations$origin)))){
    category_<-unique(Combinations$origin)[j] 
    sub<-which(Combinations$origin==category_)
    SUBD_df<-Combinations[sub,]
    colnames(SUBD_df)[1:2]<-c("X_regressor","Y_target")
    if(j==1){
      SUBD_df<-SUBD_df %>% group_by(Tau,group,origin)%>%
        arrange(X_regressor) 
      ggdeb<-ggplot2::ggplot(data=Obs_ij,aes(x=X_regressor,
                                    y=Y_target))+
        xlab(xlab_title)+ylab(ylab_title)+
        geom_point(aes(shape="data"),
                   col="#a9a9a9",size=0.8)
      GG_QREG_biv<-ggdeb+
        geom_line(data=SUBD_df,
                  aes(x=X_regressor,y=Y_target,
                      group=interaction(Tau,origin,group),
                      col=group))+
        facet_wrap(~Tau,
                   scales = "free_y")+
        geom_ribbon(data=SUBD_df,aes(ymin=ymin,ymax=ymax,
                  group=interaction(Tau,group),
                  fill=group,
                  col="confidence_band"),
                  alpha=0.2,linetype="dashed") 
    }else{
      SUBD_df<-SUBD_df %>% group_by(Tau,group,origin)%>%
        arrange(Y_target) 
      ggdeb<-ggplot2::ggplot(data=Obs_ij,aes(y=X_regressor,
                                    x=Y_target))+
        ylab(xlab_title)+xlab(ylab_title)+
        geom_point(aes(shape="data"),
                   col="#a9a9a9",size=0.8)
      GG_QREG_biv<-ggdeb+
        geom_line(data=SUBD_df,
          aes(x=Y_target,y=X_regressor,
              group=interaction(Tau,origin,group),
              col=group))+
        facet_wrap(~Tau,
                scales = "free_y")+
        geom_ribbon(data=SUBD_df,aes(y=X_regressor,
                x=Y_target, ymin=ymin,ymax=ymax,
                group=interaction(Tau,
                group),fill=group,
                col="confidence_band"),
                linetype="dashed",
                alpha=0.2) 
    }
    list_gg[[j]]<-GG_QREG_biv+
      guides(linetype="none",shape="none")+
      scale_color_manual(values=cols_ggplot)+
      scale_fill_manual(values=cols_ggplot)+
      labs(col="Legend",
           fill="Legend")
  }
  return(list_gg)
}
#' Bivariate return curve of observed and 
#' simulated time series for several return levels
#'
#' @param Obs_ij Matrix. Bivariate observed values.
#' @param Simul_ij Matrix. Bivariate values from 
#' simulated multivariate time series. 
#' @param l_name_time list. Variable and 
#' time step analysed.
#' @param cols_ggplot vector. Colors used 
#' in the ggplot used.
#' @param l_unit vector[string]. Variable unit.
#' @param levels_used Float. Level of the return curve.
#' @param zeta_sim float. Modification of the 
#' return period to consider the effect of 
#' simulating specifically extreme time series.
#' @param hearts object from the parallel package.
#' @param sub_pos Boolean. Select the indexes of 
#' positive values. 
#' @param Nb_boot int. Number of bootstrap samples 
#' to provide a confidence band.
#' @param Obs_ext_ij Matrix. Bivariate values from 
#' observed extreme multivariate time series.
#'
#' @return Ggplot object. 
#' @export
#'
#' @examples
Bivariate_RLevel_simul_vs_obs<-function(Obs_ij,Simul_ij,
                                        l_name_time,cols_ggplot,l_unit,
                                        levels_used,zeta_sim,
                                        hearts,sub_pos=TRUE,
                                        Obs_ext_ij,Nb_boot){
  
  if(sub_pos){
    ### Need to modify zeta if positive values are used.
    ### Use Bayes to retrieve link extremes/positive values
    
    d<-ncol(Obs_ij)
    Sum_obs<-rowSums(Obs_ij>0)
    Sum_cond_obs<-rowSums(Obs_ext_ij>0)
    Ind_cond_obs<-which(Sum_cond_obs==d)
    
    Sum_sim<-rowSums(Simul_ij>0)
    Ind_obs<-which(Sum_obs==d)
    
    Numerator<-length(Ind_cond_obs)
    Denominator<-length(Ind_obs)
    
    ### from Bayes formula
    zeta_sim<-Numerator/Denominator
    Ind_sim<-which(Sum_sim==d)
    ### Select
    Obs_ij<-Obs_ij[Ind_obs,]
    Simul_ij<-Simul_ij[Ind_sim,]
  }
  ### Empirical bivariate return levels
  whill <- seq(0, 1, by = 0.001)
  expdata<-ReturnCurves::margtransf(data = Obs_ij, 
                                    qmarg = rep(0.95, 2),
                                    constrainedshape = T)
  Curves_found<-lapply(levels_used,FUN = function(x){
    
    rch<-ReturnCurves::rc_est(margdata = expdata, w = whill, 
                              p = x, method = "hill",
                              q = 0.95, constrained = F)
    rch_unc<-ReturnCurves::rc_unc(rch, nboot = Nb_boot, 
                                  nangles = 200, alpha = 0.05)
    Curve<-rch_unc@retcurve@rc
    Lower<-rch_unc@unc$lower
    Upper<-rch_unc@unc$upper
    bounds_type<-c(rep("Qinf",nrow(Lower)),
                   rep("Qsup",nrow(Upper)))
    Bounds<-rbind.data.frame(Lower,
                             Upper)
    Tau_<-rep(x,nrow(Curve))
    Curve_modif<-cbind.data.frame(Curve,Tau_)
    Tau<-rep(x,nrow(Bounds))
    Df<-cbind.data.frame(Tau,Bounds,bounds_type)
    return(list("Unc"=Df,
                "Curve"=Curve_modif))
  })
  Fct_extract_curve<-function(x,name_key){
    return(x[[name_key]])
  }
  List_Curve<-lapply(Curves_found,Fct_extract_curve,
                     name_key="Curve")
  List_Unc<-lapply(Curves_found,Fct_extract_curve,
                   name_key="Unc")
  combined_df<-do.call(rbind.data.frame,  
                       List_Curve)
  Bounds_for_graph<-do.call(rbind.data.frame,  
                            List_Unc)
  colnames(combined_df)<-c("t","s","level")
  combined_df$level<-as.character(combined_df$level)
  
  # Bounds_for_graph<-do.call(rbind.data.frame,Rult)
  colnames(Bounds_for_graph)<-c("level","t","s","bounds")
  Bounds_for_graph$level<-paste0("rho==",
                                 as.character(Bounds_for_graph$level))
  poly_data <- Bounds_for_graph %>%
    group_by(level) %>%
    group_modify(~{
      
      lower <- .x %>%
        filter(bounds == "Qinf") %>%
        arrange(t)
      
      upper <- .x %>%
        filter(bounds == "Qsup") %>%
        arrange(desc(t))
      
      bind_rows(lower, upper)
      
    })
  
  ### Curves for simulated
  ### consider formula of total probs.
  levels_used_for_sim<-levels_used/zeta_sim

  Curves_found_SIM<-lapply(levels_used_for_sim,
       FUN = function(x){
         JCsim<-texmex::JointExceedanceCurve(Sample = Simul_ij,
                                             ExceedanceProb = x)
         return(JCsim)
         
       })
  
  All_sim<-lapply(c(1:length(levels_used)),function(x){
    df_j<-do.call(cbind.data.frame,Curves_found_SIM[[x]])
    nj<-nrow(df_j)
    ### Add levels used
    name_j<-rep(levels_used[x],nj)
    df_j$cl<-name_j
    return(df_j)
  })
  combined_df_sim<- do.call(rbind, All_sim)
  colnames(combined_df_sim)<-c("t","s","level")
  combined_df_sim$level<-paste0("rho==",
                                as.character(combined_df_sim$level))
  ### Curves for data Texmex
  Texmex_curves<-list()
  d<-length(l_name_time)
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
  hh<-lapply(Curves_tex_,FUN = function(x){
    level<-attributes(x)$ExceedanceProb
    Vectors_x_y<-cbind(x[[1]],x[[2]])
    LEVEL<-rep(level,nrow(Vectors_x_y))
    return(cbind(Vectors_x_y,LEVEL))
  })
  combined_df_TEX<-do.call(rbind.data.frame,hh)
  colnames(combined_df_TEX)<-c("t","s","level")
  combined_df_TEX$level<-paste0("rho==",
                                as.character(combined_df_TEX$level))
  ### Create ggplot object
  custom_labels <-as.character(levels_used)
  GG_BivRL_obs_with_Tawn<-ggplot2::ggplot(Obs_ij)+
    geom_point(aes(x=V1,y=V2,col="data"),size=0.5)
  
  exp_type<-expression(x[M]^t~v)
  NAMES<-names(l_name_time)
  N1<-NAMES[1]
  N2<-NAMES[2]
  N1_graph<-Fct_correct_name(x = N1,target = "Surcote",
                             replacement = "Surge")
  N2_graph<-Fct_correct_name(x = N2,target = "Surcote",
                             replacement = "Surge")
  T1<-l_name_time[[N1]]
  T2<-l_name_time[[N2]]
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t=T1,
                                              x=N1_graph,v=l_unit[[N1]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                                         list(t=T2,
                                              x=N2_graph,v=l_unit[[N2]])))
  combined_df$level<-paste0("rho==",
                            combined_df$level)
  GG_BivRL_obs<-GG_BivRL_obs_with_Tawn+
    geom_line(data=combined_df,aes(x=t,y=s,
                                   group=interaction(level),
                                   col="data"))+
    facet_wrap(~level,
               labeller = label_parsed)+
    geom_line(data=combined_df_TEX,aes(x=t,y=s,
                                       group=interaction(level),
                                       col="CEVmodel"))+
    geom_line(data=combined_df_sim,aes(x=t,y=s,
                                       group=interaction(level),
                                       col="simulations"))+
    geom_polygon(
      data = poly_data,
      aes(x = t, y = s,  
          group = interaction(level)),fill="lightgrey",
      alpha = 0.3
    )+
    geom_line(data = Bounds_for_graph,
              aes(x = t, y = s,
                  group = interaction(level, bounds),
                  col="confidence_band"),
              linetype="dashed"
    ) +
    labs(col="Legend")+
    guides(fill="none")+
    xlab(xlab_title)+
    ylab(ylab_title)+
    theme(axis.title=element_text(size=17),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text=element_text(size=14))+
    scale_color_manual(values = cols_ggplot)
  
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

Marg_2d_cond_ext<-function(Simul_ij,Obs_ij,l_name_time,
                               cols_ggplot,l_unit,risk_function,
                           ind_var_cond){
  NAMES<-names(l_name_time)
  ### Select knowing that one coordinate is above the threshold
  
  Th_min<-as.numeric(quantile(Obs_ij[,ind_var_cond],0.50))
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
  N1<-NAMES[1]
  N2<-NAMES[2]
  N1graph<-Fct_correct_name(x = N1,target = "Surcote",
                            replacement = "Surge")
  N2graph<-Fct_correct_name(x = N2,target = "Surcote",
                            replacement = "Surge")
  T1<-l_name_time[[N1]]
  T2<-l_name_time[[N2]]
  xlab_title<-do.call("substitute", list(exp_type[[1]], 
                           list(t=T1,
                                x=N1graph,v=l_unit[[N1]])))
  ylab_title<-do.call("substitute", list(exp_type[[1]], 
                           list(t=T2,
                                x=N2graph,v=l_unit[[N2]])))
  GG0<-ggplot2::ggplot(data=Sub_obs_ij,
            aes(x=V_1,y=V_2,col="data"))+
    geom_point(size=2,pch=19)+
    geom_point(data=Sub_sim_ij,
            aes(x=V_1,y=V_2,col="simulations"),
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




