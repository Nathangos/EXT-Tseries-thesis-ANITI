F_automatic_ggplot_marg<-function(vect_time_choice, 
                                  l_unit,risk_function,
                                  ind_var_cond){
  obj<-Marg_2d_simul_vs_obs(list_simul = L_simul,
                            list_obs = L_obs,
                            vect_times = vect_time_choice,
                            l_name_2V =L_nameV ,
                            cols_ggplot=cols_,
                            l_unit = l_unit)
  ggsave(filename = paste0(TEND_path
                           ,"/",Prefix,"_distrib_values_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,
                           opt_used,".png"),
         width=8,height = 6,
         plot = obj)
  ## Cond_ext
  obj_cond<-Marg_2d_cond_ext(list_simul = L_simul,
                             list_obs = L_obs,
                             vect_times = vect_time_choice,
                             l_name_2V =L_nameV ,
                             cols_ggplot=cols_,
                             l_unit = l_unit,
                             risk_function =risk_function,
                             ind_var_cond = ind_var_cond)
  ggsave(filename = paste0(TEND_path
                           ,"/",Prefix,"_cond_distrib_values_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,
                           opt_used,".png"),
         width=8,height = 6,
         plot = obj_cond)
}
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

Extremo_Quality_MV<-function(list_simul,list_reality,name_cond,name_other,q_per_time
                             ,B){
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
  CROSS_extremo_sim<-cross_extremogram(Matrix_pairs = Matrix_pairs,
                                       var_cond = Variable_simul_cond,
                                       var_other = Variable_simul_other,
                                       l_Tau = list_Tau)
  
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
  result_delta<-df_extremo %>% group_by(time_d) %>% summarise(val_data=mean(reality_extremo),
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
