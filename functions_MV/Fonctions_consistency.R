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
Fct_correlations<-function(method_corr,df1,
                           Intersect_times,df2=NA,
                           mat_corr=TRUE){
  if(is.null(dim(df2))){
    Na_found<-(!is.na(df2))
    Cond<-sum(Na_found)==length(df2)
  }else{
    Na_found<-(is.na(df2)) 
    Test<-sum(colSums(Na_found))
    Cond<-Test==0
  }
  if(Cond){
    vect_r<-sapply(X = Intersect_times,
                   FUN = function(j,x_1,x_2){
                     series_1<-x_1[,j]
                     
                     series_2<-x_2[,j]
                     return(cor(x = series_1,y = series_2,method = method_corr))
                   },x_1=df1,x_2=df2)
    return(vect_r)
  }else{
    d<-ncol(df1)
    NAMES<-colnames(df1)
    if(mat_corr){
      Mat_corr<-cor(df1,method = method_corr)
      return(Mat_corr)
    }else{
      Matrix_cases<-t(utils::combn(x = c(1:d),m = 2))
      list_corr<-list()
      list_pairs<-list()
      for(j in c(1:nrow(Matrix_cases))){
        Pair<-Matrix_cases[j,]
        Pair_used<-NAMES[Pair]
        Sub_df<-df1[,Pair_used]
        value_corr<-cor(x = Sub_df[,1],y = Sub_df[,2],
                      method = method_corr)
        list_corr[[j]]<-value_corr
        list_pairs[[j]]<-paste0("(",Pair_used[1],",",
                                Pair_used[2],")")
      }
      return(data.frame("value"=unlist(list_corr),
                        "pair"=unlist(list_pairs)))
    }
  }
  
}

Resamples_correlations<-function(l_extremes_indus,B,method_corr,
                                 l_namei,l_namej){
  N<-nrow(l_extremes_indus[[1]])
  Indexes_B<-sample(x = c(1:N),size = B,replace = TRUE)
  DF1<-l_extremes_indus[[l_namei]][Indexes_B,]
  DF2<-l_extremes_indus[[l_namej]][Indexes_B,]
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
  
  Sample_inds<-sample(c(1:nrow(df1)),size =nrow(df1),
                      replace = TRUE)
  Intersection_times<-intersect(colnames(df1),colnames(df2))
  Sample_df1<-as.data.frame(df1[Sample_inds,])
  Sample_df2<-as.data.frame(df2[Sample_inds,])
  colnames(Sample_df1)<-colnames(df1)
  colnames(Sample_df2)<-colnames(df2)
  Values_obtained<-Fct_correlations(method_corr = "kendall",
                   df1 =Sample_df1, df2=Sample_df2,
                   Intersect_times = Intersection_times)
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
### Work on distrib of Omega
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
#' Title
#'
#' @param Obj_data: df. Concatenated time series. 
#' @param Obj_sim: df. Simulated concatenated angles.
#' @param l_name_time: list[str:int]. List containing
#' for each selected variable the analysed time. 
#' @param cols_ggplot: vector[str]. Vector of desired cols in 
#' the ggplot plot. 
#' @return GGplot comparing the bivariate distribution of simulated angle 
#' and observed from extreme time series. 
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
           Simul_ij = Simul_ij,
           l_name_time = l_name_time,
           cols_ggplot = cols_,
           l_unit = l_unit,
           levels_used = levels_used,
           zeta_sim = zeta_forsim,
           hearts = hearts,
           Obs_ext_ij=Obs_with_cond,
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
  Seq_th<-seq.int(from = 0.5,to = 0.95,length.out =2)
  list_ggplot_qreg<-Bivariate_Qreg_simul_vs_obs(Obs_ij = Obs_without_cond,
          Simul_ij = Simul_ij,l_name_time=l_name_time,
          cols_ggplot = cols_,l_unit = l_unit,
           zeta_sim = zeta_forsim,basis_functionGAM = basis_functionGAM,
           hearts = hearts,n_sim=n_sim,Seq_th=Seq_th)
  L_gg_plot<-length(list_ggplot_qreg)
  for(j in c(1:L_gg_plot)){
    plot_elt<-list_ggplot_qreg[[j]]+
      theme_common
    F_namej<-paste0(Path_EXT,"/",Prefix_export,
                    "qreg_",N1,"t=",
                    T1,"X",N2,"t=",
                    T2,CPLMT_export,"_d=",D,
                    opt_used,"_",j,".png")
    ggsave(filename = F_namej,
           width=8,height = 6,
           plot = plot_elt)
  }
  
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
#' @param L_obs: list [dataframe].For each key, a dataframe 
#' of simulated time series. 
#' @param l_name_time: vector[str:int]. Vector of analysed times.
#'@param l_unit: list[str]. Unit_used
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
  # NSimul_for_compar<-nrow(list_simul[[1]])
  # NObs_for_compar<-nrow(L_obs[[1]])
  # L_asked<-length(l_name_time)
  # Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
  #                  ncol=L_asked )
  # Obs_ij<-matrix(NA,nrow = NObs_for_compar,
  #                ncol=L_asked)
  # Names<-names(l_name_time)
  # for(j in c(1:L_asked)){
  #   NJ<-Names[j]
  #   if(NJ%in%names(list_simul)){
  #     value_t<-l_name_time[[j]]
  #     Simul_namev<-list_simul[[NJ]][,value_t]
  #     Simul_ij[,j]<-Simul_namev
  #     Obs_namev<-L_obs[[NJ]][,value_t]
  #     Obs_ij[,j]<-Obs_namev
  #   }
  # }
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
    # texmex::JointExceedanceCurve(Sample = Obs_ij,
    #                      ExceedanceProb = x)
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
  
  # return(list("Unc"=Unc,
  #             "Curve"=Curves))
  # return(list("Curves"=List_Curve,
  #             "Unc"=List_Unc))
  # Curves_found<-cbind.data.frame(Curves_found)
  # All<-lapply(c(1:length(levels_used)),function(x){
  #   df_j<-data.frame(Curves_found[[x]])
  #   nj<-nrow(df_j)
  #   name_j<-rep(levels_used[x],nj)
  #   df_j$cl<-name_j
  #   return(df_j)
  # })
  # combined_df <- do.call(rbind, All)
  # colnames(combined_df)<-c("t","s","level")
  # combined_df$level<-as.character(combined_df$level)
  ### Bootstrap confidence regions
  # Result_boot_confRegions<-replicate(,
  #             Onesample_Bootstrap_Bivar_RL(Obs_t_s = Obs_ij,
  #                                          levels_used = levels_used))
  # Df_transf<-apply(X=Result_boot_confRegions,
  #                  as.data.frame,MARGIN=2)
  # m<-150
  # vect_ind<-c(1:m)
  # Theta_vector_target<-pi*(m+1-vect_ind)/(2*(m+1))
  # Rult<-lapply(Theta_vector_target,Candidates_per_theta,
  #              list_boot_samples=Df_transf,
  #              hearts=hearts)
  # Bounds_for_graph<-do.call(rbind.data.frame,Rult)
  colnames(Bounds_for_graph)<-c("level","t","s","bounds")
  Bounds_for_graph$level<-paste0("rho==",
                          as.character(Bounds_for_graph$level))
  # Build polygon data
  # x and y are changing so we must create a ggplot polygon object
  # arrange(level,t) for ordering
  # group_by level--> polygon per level
  # bind_rows to assemble rows per group
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
  # expsim<-ReturnCurves::margtransf(data = Simul_ij, 
  #                          qmarg = rep(0.95, 2),
  #                          constrainedshape = T)
  Curves_found_SIM<-lapply(levels_used_for_sim,
    FUN = function(x){
    JCsim<-texmex::JointExceedanceCurve(Sample = Simul_ij,
                                 ExceedanceProb = x)
    return(JCsim)
    # rch_sim<-ReturnCurves::rc_est(margdata = expsim, 
    #                 w = whill, 
    #                 p = x, method = "hill",
    #                 q = 0.95, constrained = F)
    # rch_unc_sim<-ReturnCurves::rc_unc(rch_sim, 
    #             nboot = 50, 
    #              nangles = 150, 
    #              alpha = 0.05)
    # Curve_sim<-rch_unc_sim@retcurve@rc
  })
  #Curves_found_SIM<-cbind.data.frame(Curves_found_SIM)
  #print(Curves_found_SIM)
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
  # expCEV<-ReturnCurves::margtransf(data = Simuls_from_tex_all, 
  #                        qmarg = rep(0.95, 2),
  #                        constrainedshape = T)
  Curves_tex_<-lapply(levels_used,FUN = function(x){
    JRC_j<-texmex::JointExceedanceCurve(
      Sample = Simuls_from_tex_all,
       ExceedanceProb = x,
      which=c("V1","V2"))
    # rch_CEV<-ReturnCurves::rc_est(margdata = expCEV,
    #                 w = whill,
    #                 p = x, method = "hill",
    #                 q = 0.95, constrained = F)
    # rch_unc_CEV<-ReturnCurves::rc_unc(rch_CEV,
    #             nboot = 50,
    #              nangles = 150,
    #              alpha = 0.05)
    # Curve_CEV<-rch_unc_CEV@retcurve@rc
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
  # return(list("dataObs"=combined_df,
  #             "Bounds"=Bounds_for_graph,
  #             "combi_sim"=combined_df_sim,
  #             "obs"=Obs_ij,
  #             "TEX"=combined_df_TEX)
  # )
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
Marg_2d_cond_ext<-function(Simul_ij,Obs_ij,l_name_time,
                               cols_ggplot,l_unit,risk_function,
                           ind_var_cond){
  
  # Z<-1
  # NSimul_for_compar<-nrow(list_simul[[1]])
  # NObs_for_compar<-nrow(L_obs[[1]])
  # Simul_ij<-matrix(NA,nrow = NSimul_for_compar,
  #                  ncol =length(vect_times) )
  # Obs_ij<-matrix(NA,nrow = NObs_for_compar,
  #                ncol =length(vect_times) )
  # for(j in c(1:length(list_simul))){
  #   if(names(list_simul)[j]%in%l_name_2V){
  #     Simul_namev<-list_simul[[j]][,vect_times[Z]]
  #     Simul_ij[,Z]<-Simul_namev
  #     Obs_namev<-L_obs[[j]][,vect_times[Z]]
  #     Obs_ij[,Z]<-Obs_namev
  #     Z<-Z+1
  #   }
  # }
  # Obs_ij<-as.data.frame(Obs_ij)
  # colnames(Obs_ij)<-sapply(c(1:length(l_name_2V)),function(x){
  #   return(paste0("V_",x))
  # })
  # Simul_ij<-as.data.frame(Simul_ij)
  # colnames(Simul_ij)<-colnames(Obs_ij)
  
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
Analysis_diag_HTawn_Conv_Hull<-function(Theta_opt,vect_name_variables){
  d<-length(Theta_opt)
  list_GG_for_hull<-list()
  vect_name_corrected<-sapply(vect_name_variables,FUN = Fct_correct_name,
         replacement="Surge",target="Surcote")
  for(Z in c(1:d)){
    HTawn_z<-Theta_opt[[Z]]
    Modelboot_z<-texmex::bootmex(x = HTawn_z)
    ### Boot run
    #j<-1
    Results<-lapply(Modelboot_z$boot,function(x){
      y<-x$dependence
      matrix<-as.data.frame(t(y[1:2,]))
      colnames(matrix)<-c("a","b")
      rownames(matrix)<-colnames(y)
      df <- tibble::rownames_to_column(matrix
                                  , "variable")
      df$model<-rep(as.character(Z),
                    nrow(df))
      return(df)
    })
    list_GG_for_hull[[Z]]<-do.call(rbind.data.frame,
                              Results)
  }
  Combinaison_for_hull<-do.call(rbind.data.frame,
                  list_GG_for_hull)
  Combinaison_for_hull$variable<-sapply(Combinaison_for_hull$variable,
              FUN = function(x){
                return(Fct_correct_name(x = x,target = "Surcote",
                                        replacement = "Surge"))
              })
  Vector_shown_chull<-sapply(vect_name_variables,
    FUN = function(x){
        return(Fct_correct_name(x = x,
                                target = "Surcote",
                                replacement = "Surge"))
  })
  DEFAULT_colors<-scales::hue_pal()(
      length(Vector_shown_chull))
  Real_labels <- parse(
    text = paste0("u[", Vector_shown_chull, "]")
  )
  Chull<-Combinaison_for_hull %>% 
    group_by(model,variable)%>%  
    slice(chull(a, b))
  
  GG_hull<-ggplot2::ggplot(Combinaison_for_hull,
                  aes(a, b)) +
    geom_point() +
    facet_wrap(~factor(variable,
              levels = vect_name_corrected))+
    geom_polygon(data =Chull,
                 alpha = 0.3,aes(fill=model)) +
    xlim(0,1)+ylim(-1,1)+
    theme_minimal()+labs(fill="Legend")+
    theme(axis.title=element_text(size=25),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          axis.text = element_text(size=12),
          strip.text = element_text(size = 13),
          legend.direction = "horizontal",
          legend.position = "bottom")
  GG_hull<-GG_hull+
    scale_fill_manual(values=DEFAULT_colors,
                      labels=Real_labels)
  return(GG_hull)
}



Analysis_diag_HTawn_evol_DQU<-function(vect_dqu,
                                       Vect_obs,MQU,
                                       chosen_dqu){
  LNAME<-colnames(Vect_obs)
  LNAME_corrected<-sapply(LNAME,
      FUN =function(x){
        Fct_correct_name(x = x,target = "Surcote",
                         replacement = "Surge")})
  d<-ncol(Vect_obs)
  GG_theta<-ggplot2::ggplot()
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
    Rult2$model<-rep(cond_var,
                     nrow(Rult2))
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
    Melting_indep$variable_cond<-rep(cond_var,
                           nrow(Melting_indep))
    Melting_indep$dqu<-rep(chosen_dquj,
                           nrow(Melting_indep))
    Indep_whole<-rbind(Indep_whole,
                       Melting_indep)
  }
  Dmodel<-length(unique(Rult_simplified$model))
  Real_test<-expression(u[j])
  if(Dmodel==1){
    Real_label<-as.expression(do.call("substitute", 
                  list(Real_test[[1]], 
                        list(j=1))))
    Rult_simplified$model<-as.character(Rult_simplified$model)
    Rult_simplified$variable<-sapply(Rult_simplified$variable,
                  FUN =function(x){
                        Fct_correct_name(x = x,target = "Surcote",
                            replacement = "Surge")})
    GG_theta<-ggplot2::ggplot(Rult_simplified
                     ,aes(x=quantile,y=value,
                          group=interaction(variable,param,model),
                          col=param))+
      geom_line()+
      ylab("Estimator")+
      xlab("Dqu")+
      facet_wrap(~variable,
                 scales = "free_y")+
      #geom_point(aes(shape=model))+
      labs(col="Legend",
           linetype="Model")+
      geom_vline(aes(xintercept=dqu,
                     linetype="dep_quantile"))+
      scale_linetype_manual(values=2,
                            labels=Real_label)
    
    Indep_whole$variable_cond<-as.character(Indep_whole$variable_cond)
    
    Indep_whole$variable<-sapply(Indep_whole$variable,
           FUN =function(x){
             Fct_correct_name(x = x,target = "Surcote",
                              replacement = "Surge")})
    GG_Indep<-ggplot2::ggplot(Indep_whole,
                     aes(x=quantile,y=pval,
                         group=interaction(variable_cond,
                                           test_used,variable),
                         col=test_used))+
      geom_line()+
      facet_wrap(~variable)+
      #geom_point(aes(shape=variable_cond))+
      ylab("p value")+
      xlab("Dqu")+
      labs(col="Legend",
           linetype="Model")+
      geom_hline(yintercept = 0.05,col="black")+
      geom_vline(aes(xintercept = dqu,
                     linetype="dep_quantile"))+
      scale_linetype_manual(values=2,
                            labels=Real_label)
  }else{
    GG_theta<-list()
    GG_Indep<-GG_theta
    for(J in c(1:Dmodel)){
      Real_label<-as.expression(do.call("substitute", 
              list(Real_test[[1]], 
                list(j=LNAME_corrected[J]))))
      Rult_simplified_i<-Rult_simplified[which(Rult_simplified$model==J),]
      Rult_simplified_i$model<-as.character(Rult_simplified_i$model)
      Rult_simplified_i$variable<-as.character(Rult_simplified_i$variable)
      Rult_simplified_i$variable<-sapply(Rult_simplified_i$variable,
             FUN =function(x){
               Fct_correct_name(x = x,target = "Surcote",
                                replacement = "Surge")})
      
      GG_theta_J<-ggplot2::ggplot(Rult_simplified_i
                       ,aes(x=quantile,y=value,
                            group=interaction(variable,param,model),
                            col=param))+
        geom_line()+
        ylab("Estimator")+
        xlab("Dqu")+
        facet_wrap(~variable,
                   scales = "free_y")+
        #geom_point(aes(shape=model))+
        labs(col="Legend","linetype"="Model")+
        geom_vline(aes(xintercept=dqu,
                       linetype="dep_quantile"))+
        scale_linetype_manual(values=2,
                              labels=Real_label)+
        theme(axis.title = element_text(size=15),
              strip.text=element_text(size=13),
              legend.title = element_text(size=14),
              legend.text=element_text(size=13),
              legend.direction = "horizontal",
              legend.position = "bottom")+
        ggguides::legend_order_guides(linetype = 1, colour = 2)
      GG_theta[[J]]<-GG_theta_J
      
      Indep_whole_i<-Indep_whole[which(Indep_whole$variable_cond==J),]
      Indep_whole_i$variable_cond<-as.character(Indep_whole_i$variable_cond)
      Indep_whole_i$variable<-as.character(Indep_whole_i$variable)
      Indep_whole_i$variable<-sapply(Indep_whole_i$variable,
         FUN =function(x){
           Fct_correct_name(x = x,target = "Surcote",
                            replacement = "Surge")})
      GG_Indep_J<-ggplot2::ggplot(Indep_whole_i,
                       aes(x=quantile,y=pval,
             group=interaction(variable_cond,
                               test_used,variable),
             col=test_used))+
        geom_line()+
        facet_wrap(~variable)+
        #geom_point(aes(shape=variable_cond))+
        ylab("p value")+
        xlab("Dqu")+
        labs(col="Legend","linetype"="Model")+
        geom_hline(yintercept = 0.05,col="black")+
        geom_vline(aes(xintercept = dqu,
                       linetype="dep_quantile"))+
        scale_linetype_manual(values=2,
                              labels=Real_label)+
        theme(legend.direction = "horizontal",
              legend.position = "bottom",
              axis.title = element_text(size=15),
              strip.text=element_text(size=13),
              legend.title = element_text(size=14),
              legend.text=element_text(size=13))+
        ggguides::legend_order_guides(linetype = 1, colour = 2)
      GG_Indep[[J]]<-GG_Indep_J
    }
  }

  return(list("theta"=GG_theta,
              "indep"=GG_Indep))
}
