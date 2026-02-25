
Create_list_fromALL<-function(l_name,obj_source){
  list_toreturn<-list()
  for(j in c(1:length(l_name))){
    if(typeof(obj_source)=="list"){
      list_toreturn[[l_name[j]]]<-obj_source[[j]]
    }
    if((typeof(obj_source)=="double")&&
       (!is.null(dim(obj_source)))){
      list_toreturn[[l_name[j]]]<-obj_source[,j]
    }
    if((typeof(obj_source)=="double")&&
       (is.null(dim(obj_source)))){
      list_toreturn[[l_name[j]]]<-obj_source[j]
    }
    
  }
  return(list_toreturn)
} 
Choice_automatic_thresh_per_variable_time<-function(var_t,nb_threshs){
  Qtiles_candidates<-seq.int(0,0.95,
                             length.out=nb_threshs)
  thresholds1 <- quantile(var_t,Qtiles_candidates)
  Metric_result<-eqd(var_t, 
                     thresh = thresholds1)
  plot(thresholds1, Metric_result$dists, 
       xlab="Threshold", ylab="Metric value")
  index_used<-which(Metric_result$thresh==thresholds1)
  return(Qtiles_candidates[index_used])
}
Choice_automatic_thresh_per_variable<-function(df,
                                               nb_threshs){
  return(apply(X = df,MARGIN = 2,FUN =
                 Choice_automatic_thresh_per_variable_time,
               nb_threshs=nb_threshs))
}

#' MarTransfo_TS_exts_Mixture
#'
#' @param type_donnees 
#' @param coeurs 
#' @param lien_racine 
#' @param liste_noms 
#' @param file_dates 
#' @param n.dens 
#' @param opt_Frech 
#' @param type_entree 
#'
#' @return
#' @export
#'
#' @examples
MarTransfo_TS_exts_Mixture<-function(type_donnees,coeurs,lien_racine,liste_noms,
                                     file_dates,n.dens,
                                     opt_Frech,type_entree,
                                     Nb_Threshs){
  p_U<-list()
  l_ALL<-list()
  l_Orig<-list()
  j<-1
  l_AD<-list()
  for(name_variable in liste_noms){
    lien_donnees<-paste0(lien_racine,name_variable,"_residuals.csv")
    if(type_entree=="clust"){
      All<-read.csv(file=lien_donnees)[,2:39]
      colnames(All)<-c(1:38)
      Vrais_indices<-All[,38]
      Donnes<-All[,c(1:37)]
    }
    else{
      Donnes<-read.csv(file=lien_donnees)[,2:38]
    }
    print(nrow(Donnes))
    colnames(Donnes)<-c(1:37)
    rownames(Donnes)<-c(1:nrow(Donnes))
    l_Orig[[name_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    ### Choice of threshold using metric
    #######
    p_U[[name_variable]]<-Choice_automatic_thresh_per_variable(df = Donnes,
                                                               nb_threshs = Nb_Threshs)
    # AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,
    #                           fonction_seuil = p_U[[name_variable]],n.dens = n.dens,
    #                           nom=name_variable,type_entree=type_entree,
    #                           dates_prises=dates_import,j_show = 19)
    # l_AD[[name_variable]]<-AN_GPD
    if(type_entree=="donnees_brutes"){
      return(TRUE)
    }
    Vecteur_DIM<-c(1:ncol(Donnes))
    Resultat_P<-as.data.frame(t(sapply(Vecteur_DIM,f_marginales_all_Pareto,data_to_tf=Donnes,
                                       p_u=p_U[[name_variable]],
                                       n.dens=n.dens)))
    K<-Resultat_P$kernel_dens
    PARETO<-Resultat_P$obs
    L_EVT<-list("gamma"=Resultat_P$gamma,
                "scale"=Resultat_P$scale,
                "threshold"=Resultat_P$threshold,
                "p_u"=Resultat_P$p_u)
    Vtransf<-do.call(cbind.data.frame,PARETO)
    if(opt_Frech==TRUE){
      Vunif<-(1-(1/Vtransf))
      Vtransf<--1/log(Vunif)
    }
    l_ALL[[name_variable]]<-list("transf"=Vtransf,"K"=K,
                                 "LEVT"=L_EVT)
    if(name_variable==liste_noms[1]){
      df<-matrix(NA,ncol=length(liste_noms),
                 nrow=nrow(Vtransf))
    }
    df[,j]<-apply(X = Vtransf,MARGIN = 1,FUN = calcul_norme_L2)
    j<-j+1
  }
  # Pareto ------------------------------------------------------------------
  return(list("resume"=l_ALL,"df"=df,"GPD"=l_AD,
              "orig"=l_Orig,"type_transfo"="Mixture_emp_GPD"))
}
fct_extract_Transf<-function(time_list){
  return(time_list[["unif_convert_t"]])
}
fct_extract_ParamsEGPD<-function(time_list){
  return(time_list[["params_egpd"]])
}
fct_extract_InitEGPD<-function(time_list){
  return(time_list[["Init"]])
}
fct_extract_FittingEGPD<-function(time_list){
  return(time_list[["Model_t"]])
}
#' Convert_time_z
#'
#' @param z: int. Time index. 
#' @param Data_pos: dataframe. Time series of positive marginals. 
#' @param show_EGPD: Bool. Show or not the details of EGPD fitting. 
#'
#' @return
#' @export
#'
#' @examples
Convert_time_z<-function(z,Data_pos,show_EGPD,
                         list_params_EGPD){
  print(paste0("Result for time ",z))
  n.cyc<- list_params_EGPD[["n.cyc"]]
  mu.step<- list_params_EGPD[["mu.step"]]
  sigma.step<- list_params_EGPD[["sigma.step"]]
  nu.step<- list_params_EGPD[["nu.step"]]
  tau.step<- list_params_EGPD[["tau.step"]]
  EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")
  Col_t<-Data_pos[,z]
  Col_pos_t<-Col_t[which(Col_t>0)]
  db<-as.data.frame(Col_pos_t)
  colnames(db)<-c("x")
  con <- gamlss.control(n.cyc = n.cyc,mu.step = mu.step,
                        sigma.step = sigma.step, nu.step =nu.step,
                        tau.step = tau.step,autostep=TRUE,
                        trace = show_EGPD)
  con.i<-glim.control(glm.trace = FALSE)
  # Identity to enable mu/gamma<0
  Th<-quantile(x =Col_pos_t,0.90)
  #Th<-0
  #FIT_pos<-mev::gp.fit(xdat = Col_pos_t,threshold =Th)$est
  FIT_pos<-extRemes::fevd(x = Col_pos_t,threshold = Th,
                          method = "Lmoments",
                          type = "GP")$results
  # GG_shapeorig<-Graphics_estimators_gamma(series = Col_pos_t,
  #                           vect_k = c(50:300),
  #                           Title_graphic = paste0("result time ",z))
  # print(GG_shapeorig)
  
  Shape<-FIT_pos[2]
  # ## Use gpd property to get the scale at 0.
  if(Shape<0){
    Scale<-FIT_pos[1]-Th*Shape
  }else{
    Scale<-FIT_pos[1]
  }
  INIT<-c(Shape,Scale)
  #INIT<-c(FIT_pos[2],FIT_pos[1])
  print(INIT)
  #INIT<-c(FIT_pos[2],FIT_pos[1])
  ## Nu.start initialisation using moments.
  Th_beg<-quantile(x = Col_pos_t,0.05)
  Lower_tail<-Col_pos_t[which(Col_pos_t<Th_beg)]
  Moment_1<-mean(Lower_tail)
  nu.start<-as.numeric((1-(Moment_1/Th_beg))^(-1)-1)
  Fitting_time_t <- gamlss(x~1, 
                           data=db,family = EGPD1Family(mu.link = "identity"),
                           control = con,mu.start=INIT[1],
                           sigma.start=INIT[2],nu.start=nu.start,
                           i.control=con.i,method=CG())
  muFit <-fitted(Fitting_time_t,"mu")[1]
  sigmaFit <- predict(Fitting_time_t,what="sigma", 
                      type="response")[[1]]
  nuFit <- predict(Fitting_time_t,what="nu", 
                   type="response")[[1]]
  Params_egpd<-list("mu"=muFit,"sigma"=sigmaFit,"nu"=nuFit)
  vect_unif_conver_t<-pEGPDModel1(q = Col_t,mu =muFit,
                                  sigma = sigmaFit,
                                  nu = nuFit)
  return(list("params_egpd"=Params_egpd,
              "unif_convert_t"=vect_unif_conver_t,
              "Model_t"=Fitting_time_t,
              "Init"=c(INIT,nu.start)))
}

MarTransfo_TS_EXTGP_MV<-function(type_donnees,coeurs,lien_racine,
                                 liste_noms,opt_Frech,n.dens,type_entree,p_U,
                                 file_dates,show_EGPD,
                                 list_params_EGPD){
  l_ALL<-list()
  l_Orig_all<-list()
  l_Orig<-list()
  l_AD<-list()
  LIST_SUB_orig<-list()
  LIST_SUB_TRANSFO<-list()
  j<-1
  df<-list()
  for(name_variable in liste_noms){
    lien_donnees<-paste0(lien_racine,name_variable,"_residuals.csv")
    if(type_entree=="clust"){
      All<-read.csv(file=lien_donnees)[,2:39]
      colnames(All)<-c(1:38)
      Vrais_indices<-All[,38]
      Donnes<-All[,c(1:37)]
    }
    else{
      Donnes<-read.csv(file=lien_donnees)[,2:38]
    }
    colnames(Donnes)<-c(1:37)
    rownames(Donnes)<-c(1:nrow(Donnes))
    l_Orig_all[[name_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    dates_import<-dates_import[-1]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,fonction_seuil = p_U[[name_variable]],n.dens = n.dens,
                              nom=name_variable,type_entree=type_entree,
                              dates_prises=dates_import,j_show = 19)
    l_AD[[name_variable]]<-AN_GPD
    if(type_entree=="donnees_brutes"){
      return(TRUE)
    }
    
    SUM_pos<-apply(Donnes>0,MARGIN = 1,FUN = sum)
    IND_posj<-which(SUM_pos==ncol(Donnes))
    if(j==1){
      IND_select<-IND_posj
    }
    IND_select<-intersect(IND_select,IND_posj)
    j<-j+1
  }
  for(name_variable in liste_noms){
    Donnes<-l_Orig_all[[name_variable]]
    colnames(Donnes)<-c(1:37)
    rownames(Donnes)<-c(1:nrow(Donnes))
    Vecteur_DIM<-c(1:ncol(Donnes))
    # 1) Use GAM to analyse short-tail distribs ----------------------------------
    Donnes_pos<-Donnes[IND_select,]
    l_Orig[[name_variable]]<-Donnes_pos
    Unif_EGPD<-matrix(NA,nrow = length(IND_select),
                      ncol=ncol(Donnes))
    ALL_results<-lapply(c(1:ncol(Donnes)),FUN =Convert_time_z,
                        Data_pos=Donnes_pos,
                        show_EGPD=show_EGPD,
                        list_params_EGPD=list_params_EGPD)
    Unif_EGPD<-cbind.data.frame(lapply(X = ALL_results,
                                       FUN = fct_extract_Transf))
    Vect_egpd<-lapply(X = ALL_results,
                      FUN = fct_extract_ParamsEGPD)
    Vect_Init<-lapply(X = ALL_results,
                      FUN = fct_extract_InitEGPD)
    Vect_fitting<-lapply(X = ALL_results,
                         FUN = fct_extract_FittingEGPD)
    
    if(opt_Frech){
      FINAL_transfo<--1/log(Unif_EGPD)
    }else{
      FINAL_transfo<-(1/(1-Unif_EGPD))-1
    }
    l_ALL[[name_variable]]<-list("transf"=FINAL_transfo,
                                 "params_transfo"=Vect_egpd,
                                 "Fitting"=Vect_fitting,
                                 "INIT"=Vect_Init)
    df[[name_variable]]<-apply(X = FINAL_transfo,MARGIN = 1,
                               FUN = calcul_norme_L2)
    
  }
  return(list("resume"=l_ALL,"df"=as.data.frame(df),"GPD"=l_AD,
              "orig"=l_Orig,"orig_all"=l_Orig_all,
              "SUB_set"=IND_select,"type_transfo"="EXTGPD"
  ))
}


#' AD_KS_B_choice
#'
#' @param NPY_per_block : Int. Length of each block. 
#' @param Vect_l_function : vector(float). Vector of l() values for each variable
#'
#' @return Dataframe(float). P values of Frechet law for each stat test (KS, AD, CVM)
#' and variable
#' @export
#'
#' @examples
AD_KS_B_choice<-function(NPY_per_block,Vect_l_function,l_name){
  
  NRows<-nrow(Vect_l_function)
  M_max<-round(NRows/NPY_per_block)
  BLOCKS<-as.character(round(c(1:NRows)/NPY_per_block))
  Calul_max<-cbind.data.frame(Vect_l_function,BLOCKS)
  colnames(Calul_max)<-c(l_name,"block")
  Vect_input<-Calul_max %>% group_by(block) %>% 
    summarise(across(l_name, ~ max(.x, na.rm = TRUE)))
  #across to consider every name and max(.x) to realise 
  #the operation, kind of sapply
  #normalisation
  Vect_input<-as.data.frame(Vect_input[,l_name])/NPY_per_block
  delta_moment<-0.5
  
  # l law -------------------------------------------------------------------
  value_found<-colMeans(Vect_input^delta_moment)
  value_std_frechet<-factorial(-delta_moment)
  scale_frechet_d<-(value_found/value_std_frechet)^(1/delta_moment)
  scale_frechet_d
  Result_std<-matrix(NA,nrow = nrow(Vect_input),
                     ncol=ncol(Vect_input))
  Output<-rep(ncol(Vect_input)*3)
  Fill<-c()
  for(j in c(1:ncol(Vect_input))){
    Result_std[,j]<-Vect_input[,j]/scale_frechet_d[j]
    pval<-goftest::ad.test(Vect_input[,j],null = extRemes::"pevd",
                           scale=scale_frechet_d[j],shape=1,
                           loc=scale_frechet_d[j],type="GEV")$p.value
    pval2<-ks.test(Vect_input[,j],extRemes::"pevd",
                   scale=scale_frechet_d[j],shape=1,
                   loc=scale_frechet_d[j],type="GEV")$p.value
    pval3<-goftest::cvm.test(Vect_input[,j],extRemes::"pevd",
                             scale=scale_frechet_d[j],shape=1,
                             loc=scale_frechet_d[j],type="GEV")$p.value
    Fill<-c(Fill,c(pval,pval2,pval3))
  }
  CHAR<-rep(c("KS","AD","CVM"),ncol(Vect_input))
  NPY_BLOCK<-rep(NPY_per_block,ncol(Vect_input))
  DF<-cbind.data.frame(Fill,CHAR)
  colnames(DF)<-c("value","name_test")
  DF$NPY_BLOCK<-rep(NPY_per_block,nrow(DF))
  D<-ncol(Vect_input)
  Index_col<-c(sapply(X = c(1:D),
                      function(x){rep(x,nrow(DF)/D)}))
  DF$index_col<-Index_col
  return(DF)
  
}
Compar_MomentF_B<-function(NPY_per_block,Vect_l_function,l_name){
  
  NRows<-nrow(Vect_l_function)
  M_max<-round(NRows/NPY_per_block)
  BLOCKS<-as.character(round(c(1:NRows)/NPY_per_block))
  Calul_max<-cbind.data.frame(Vect_l_function,BLOCKS)
  colnames(Calul_max)<-c(l_name,"block")
  Vect_input<-Calul_max %>% group_by(block) %>% 
    summarise(across(l_name, ~ max(.x, na.rm = TRUE)))
  #across to consider every name and max(.x) to realise 
  #the operation, kind of sapply
  #normalisation
  Vect_input<-as.data.frame(Vect_input[,l_name])/NPY_per_block
  delta_moment<-0.5
  
  # l law -------------------------------------------------------------------
  value_found<-colMeans(Vect_input^delta_moment)
  value_std_frechet<-factorial(-delta_moment)
  scale_frechet_d<-(value_found/value_std_frechet)^(1/delta_moment)
  scale_frechet_d
  Result_std<-matrix(NA,nrow = nrow(Vect_input),
                     ncol=ncol(Vect_input))
  Output<-rep(ncol(Vect_input)*3)
  Fill<-c()
  for(j in c(1:ncol(Vect_input))){
    Result_std[,j]<-Vect_input[,j]/scale_frechet_d[j]
    pval<-goftest::ad.test(Vect_input[,j],null = extRemes::"pevd",
                           scale=scale_frechet_d[j],shape=1,
                           loc=scale_frechet_d[j],type="GEV")$p.value
    pval2<-ks.test(Vect_input[,j],extRemes::"pevd",
                   scale=scale_frechet_d[j],shape=1,
                   loc=scale_frechet_d[j],type="GEV")$p.value
    pval3<-goftest::cvm.test(Vect_input[,j],extRemes::"pevd",
                             scale=scale_frechet_d[j],shape=1,
                             loc=scale_frechet_d[j],type="GEV")$p.value
    Fill<-c(Fill,c(pval,pval2,pval3))
  }
  CHAR<-rep(c("KS","AD","CVM"),ncol(Vect_input))
  NPY_BLOCK<-rep(NPY_per_block,ncol(Vect_input))
  DF<-cbind.data.frame(Fill,CHAR)
  colnames(DF)<-c("value","name_test")
  DF$NPY_BLOCK<-rep(NPY_per_block,nrow(DF))
  D<-ncol(Vect_input)
  Index_col<-c(sapply(X = c(1:D),
                      function(x){rep(x,nrow(DF)/D)}))
  DF$index_col<-Index_col
  return(DF)
  
}
#' Title
#'
#' @param result_transformation : result of previous function.
#' @param l_variables : list[str]. Liste of variables.
#' @param Q_thresh : float. Proportion of extremes.
#' @param NbScores_Omega : int. Number of scores of the PCA. 
#' @param M : int. Number of simulated time series. 
#'
#' @return
#' @export
#'
#' @examples
Simul_MV_residuals<-function(result_transformation,l_variables,Q_thresh,
                             NbScores_Omega,M,rotations_available,
                             Params_risk_Function,root_for_export,
                             list_nb_scores,One_PCA_base,f_transf,
                             f_transf_inv,opt_Frech,
                             cols_ggplot){
  Shape_evd<-Params_risk_Function[["Shape_parameter"]]
  # CF Kokozka ---------------------------------------------------------------
  L<-300
  vect_k<-c(20:L)
  Array_simul<-array(data = NA,dim = c(length(l_variables),M,37))
  
  LISTE_obs_exts<-list()
  Liste_all<-result_transformation$resume
  Vect_l_function<-result_transformation$df
  # RiskF -------------------------------------------------------------------
  ### ----------------------------------------------------------------------
  for(ztilde in c(1:ncol(Vect_l_function))){
    Graphics_estimators_gamma(series = Vect_l_function[,ztilde],
                              vect_k = vect_k,
                              NB_years = 37,
                              Title_graphic = paste0("Shape parameter of l(T("
                                                     ,l_name[ztilde]
                                                     ,")"))
  }
  N<-nrow(Vect_l_function)
  L_Transf_vect_l<-list()
  L_Convert<-list()
  dir_AIC<-"RisKFunctions/evol_AIC_param"
  for(nameV in l_variables){
    #log
    VECT_lj<-log(Vect_l_function[,nameV]+1)
    # Opitz way
    # L_Transf_vect_l[[nameV]]<-f_marginales_Pareto(variable =VECT_lj,p_u = 0.05,
    #                              n.dens = 10^(4))
    # L_Convert[[nameV]]<-L_Transf_vect_l[[nameV]][["obs"]]
    # EGPD way
    Values_gp<-mev::gp.fit(VECT_lj,threshold = 0)$est
    Th_beg<-quantile(x = VECT_lj,0.10)
    Lower_tail<-VECT_lj[which(VECT_lj<Th_beg)]
    Moment_1<-mean(Lower_tail)
    nu.start<-as.numeric((1-(Moment_1/Th_beg))^(-1)-1)
    PWM_EXTGP<-mev::fit.extgp(data = VECT_lj,model = 1,init=c(nu.start,Values_gp[1],0.2),
                              method="pwm",R=20)
    Theta<-PWM_EXTGP$fit$pwm
    Name_file_plot<-paste0(root_for_export,"RisKfunctions/",nameV,"_EGPD_fitting_RiskF.png")
    MLE_EXTGP<-mev::fit.extgp(data = VECT_lj,model = 1,init=Theta,
                              method="mle",R=20)$fit$mle
    UNIF_extgp<-sapply(VECT_lj,mev::pextgp,kappa=MLE_EXTGP[["kappa"]],
                       sigma = MLE_EXTGP[["sigma"]],
                       xi=MLE_EXTGP[["xi"]])
    Quantile_LEVELS<-c(1:length(VECT_lj))/(length(VECT_lj)+1)
    ### Confidence bound
    Result_boot<-Bootstrap_conf_band(M = 500,
                                     n = length(VECT_lj),
                                     theta_egpd = MLE_EXTGP,
                                     order_quantiles = Quantile_LEVELS,
                                     alpha = 0.05)
    DF_egpd_emp<-Result_boot
    DF_egpd_emp$emp<-sort(VECT_lj)
    DF_egpd_emp$Levels<-Quantile_LEVELS
    DF_egpd_emp<-as.data.frame(DF_egpd_emp)
    GG_fit_vect_R<-ggplot(data=DF_egpd_emp,aes(x=emp,y=mean,
                                               col="model"))+
      geom_point()+
      geom_line(aes(x=emp,y=emp,col="data"))+
      geom_ribbon(aes(ymin=bound_inf,ymax=bound_sup,
                      col="confidence_band"),
                  linetype="dashed",alpha=0.15,
                  fill="grey")+
      theme(axis.title=element_text(size=20),
            legend.title = element_text(size=14),
            legend.text=element_text(size=13),
            axis.text=element_text(size=13))+
      scale_color_manual(values=cols_)+
      labs(col="Legend")+
      xlab("Empirical quantiles")+
      ylab("Fitted quantiles")
    ggsave(filename = Name_file_plot,
           plot = GG_fit_vect_R,
           width = 8,height=6)
    
    L_Transf_vect_l[[nameV]]<-MLE_EXTGP
    Frechet_j<-(-log(UNIF_extgp))**(-1)
    print(paste0("Frechet test for conversion of l(",nameV,")"))
    print(goftest::ad.test(Frechet_j,null = extRemes::"pevd",
                           shape=1,scale=1,loc=1,type="GEV")$p.value)
    L_Convert[[nameV]]<-Frechet_j
  }
  Vect_l_transf<-as.data.frame(L_Convert)
  ### Regular variations hypothesis
  Simple_MRV_HRV_analysis(couple_UV =  Vect_l_transf,
                          root_graphics = paste0(root_for_export,
                                                 "RisKfunctions/"),
                          Vect_k = c(20:500),
                          l_name = l_variables,q = 0.80)
  Vect_lunif<-exp(-Vect_l_transf^(-1))
  ### Asymmetry test 
  test_asy<-copula::exchTest(Vect_lunif,N = 1000)
  print(test_asy)
  ### What would be the copula
  Select_cop_fam<-VineCopula::BiCopSelect(Vect_lunif[,1],
                                          Vect_lunif[,2])
  FAMILY_set<-c(1:40,104,114,
                124,134,204,214,
                224,234)
  # Create a results data frame
  results_all<-sapply(FAMILY_set,FUN = FIT_cop
                      ,u1 = Vect_lunif[,1],
                      u2 = Vect_lunif[,2])
  Val<-sapply(results_all,FUN=function(x){
    return(length(x))
  })
  Ind_pos<-which(Val>0)
  List_pos<-lapply(Ind_pos,FUN=function(x){
    results_all[[x]]})
  DF<-t(cbind(sapply(List_pos,unlist)))
  DF<-as.data.frame(DF)
  DF$AIC<-round(as.numeric(DF$AIC),2)
  DF$LogLik<-round(as.numeric(DF$LogLik),2)
  DF$BIC<-round(as.numeric(DF$BIC),2)
  DF$par<-round(as.numeric(DF$par),2)
  DF$par2<-round(as.numeric(DF$par2),2)
  Sort_df<-head(DF[order(DF$AIC),c(1:6) ],5)
  write.csv(x = Sort_df,file = paste0(root_for_export,
                                      dir_AIC,"_",
                                      Approach,".csv"))
  ### 
  Vect_Pareto_margins<-apply(Vect_lunif,
                             MARGIN=2,FUN = evd::qgpd,loc=0,
                             scale=1, shape=1)
  l_prime<-apply(Vect_l_transf,MARGIN = 1,
                 FUN = List_Params_RF[["function"]])
  TS_param<-0.30
  Result_k0_metric2<-tea::mindist(data = l_prime,
                                  method = "mad",ts =TS_param )
  kstar2<-Result_k0_metric2$k0
  print(paste0("We would take ",kstar2," extreme MV time series with MAD distance"))
  Q_alt<-1-(kstar2/nrow(Vect_l_transf))
  print("Alternative proportion of extreme observations")
  print(Q_alt)
  vect_k<-c(50:200)
  GG_gam<-Graphics_estimators_gamma(series = l_prime,
                                    vect_k = vect_k,
                                    Title_graphic = "Shape parameter of gol",
                                    NB_years = 37)
  print(GG_gam)
  l_Pareto<-apply(Vect_Pareto_margins,MARGIN = 1,
                  FUN = List_Params_RF[["function"]])
  
  # Choice of the individuals -----------------------------------------------
  Seuil_lprime<-quantile(l_prime,
                         Q_thresh)
  Indices_exts<-which(l_prime>Seuil_lprime)
  Seuil_Margins_Pareto<-quantile(l_Pareto,
                                 Q_thresh)
  Excedents_lprime<-Vect_l_transf[Indices_exts,]
  
  # Export radial extreme elts ----------------------------------------------
  Ext_cat<-as.character(as.numeric(l_prime>Seuil_lprime))
  Df_rad_ext<-cbind.data.frame(Vect_l_function,Ext_cat)
  NC<-length(colnames(Vect_l_function))
  Cols_first<-sapply(c(1:NC),function(x){return(paste0("V",x))})
  colnames(Df_rad_ext)<-c(Cols_first,"category")
  GG_rad_ext_bool<-ggplot(data=Df_rad_ext,aes(x=V1,y=V2,col=category))+
    geom_point()+
    guides(col="none")+
    scale_x_continuous(transform="log10")+
    scale_y_continuous(transform="log10")+
    geom_xsidedensity(data=Df_rad_ext,aes(fill=category), alpha = 0.5)+
    geom_ysidedensity(data=Df_rad_ext,aes(fill=category), alpha = 0.5)+
    labs(fill="Legend")+
    xlab(paste0("T(",l_variables[1],") (.)"))+
    ylab(paste0("T(",l_variables[2],") (.)"))
  ggsave(filename = paste0(root_for_export,
                           "repar_points_Rad_",
                           l_variables[1],"_",
                           l_variables[2],".png"),
         plot = GG_rad_ext_bool,
         width = 8,height=6)
  
  Theta_obtained<-atan(Vect_Pareto_margins[,2]/Vect_Pareto_margins[,1])
  png(filename= paste0(root_for_export,
                       "/RisKfunctions/ang_dens_R_",
                       l_variables[1],"_",
                       l_variables[2],".png"),
      width=600,height=800)
  plot(density(Theta_obtained),
       main="Angular density of the risk functions",
       cex.lab = 1.5,
       cex.main=1.5)
  abline(v=0,col="red")
  abline(v=pi/2,col="red")
  dev.off()
  ### Distribution of angular comp for exceedances
  Exced_PTO_marg<-Vect_Pareto_margins[Indices_exts,]
  l_exts<-l_Pareto[Indices_exts]
  Angle_exceeds<-t(t(Exced_PTO_marg)%*%diag(l_exts^(-1)))
  
  df_angle_exceeds<-as.data.frame(Angle_exceeds)
  colnames(df_angle_exceeds)<-c("V1","V2")
  GG_dens_Ang<-ggplot(data = df_angle_exceeds,aes(x=V1,y=V2))+
    geom_point(col="blue")+
    geom_xsidedensity(alpha=0.5)+
    geom_ysidedensity(alpha=0.5)+
    xlab(paste0("Angular component for risk of ",
                l_name[1]))+
    ylab(paste0("Angular component for risk of ",
                l_name[2]))
  ggsave(filename = paste0(root_for_export,
                           "/RisKfunctions/Angular_val_exts_",
                           l_variables[1],"_",
                           l_variables[2],".png"),
         plot = GG_dens_Ang,width=6,height=8)
  
  Excesses_maximum<-l_prime[Indices_exts]
  # Anderson-Darling test for the maximum -- Pareto ?  ----------------------
  print("p value AD test Pareto")
  print(goftest::ad.test(x =Excesses_maximum, 
                         null = extRemes::"pevd",
                         shape=1,scale=Seuil_lprime,
                         threshold=Seuil_lprime,
                         type="GP")$p.value)
  Theta_opt<-NA
  # else{
  #   ### Drop the AD hypothesis
  #   
  #   d <- ncol(Vect_l_transf)
  #   start <- starting_point(Vect_l_transf,d)
  #   start <- c(t(start))
  #   n<-nrow(Vect_l_transf)
  #   p <- 0.3
  #   ### Number of assumed clusters.
  #   r<-2
  #   lambda <- 10^(-3)
  #   k<-(1-Q_thresh)
  #   Vect_Unit_Pareto_margins<-Vect_Pareto_margins+1
  #   if(Params_risk_Function[["name_model"]]=="HR"){
  #     points <- c(0 , 1/6,  1/8 ,  1/4 ,  1/3 , 1/2 , 2/3 , 3/4 , 1)
  #     Grid_points<- selectGrid(points, d = d, nonzero = d-1)
  #     Mod<-"SSR_row_HR"
  #   }else{
  #     ###Log variant
  #     points_log <- c(0,1/4 , 1/3 , 1/2 , 3/4 ,1)
  #     Grid_points<- tailDepFun::selectGrid(cst = points_log, d = d, 
  #                                          nonzero =d-1)
  #     Mod<-"SSR_row_log"
  #   }
  #   q<-nrow(Grid_points)
  #   w_total <- sapply(1:q, function(m) 
  #     stdfEmp(Vect_Unit_Pareto_margins, 
  #             k * n, Grid_points[m, ]))
  #   A_Theta<-param_estim(d , r , grid = Grid_points , lambda=lambda , 
  #                        num_col = NULL ,
  #                        start , type =Mod , p = p ,  
  #                        w = w_total )
  #   Vect_lambda<-c(1e-04,1e-03,1e-02,1e-01)
  #   L_Lambda<-list()
  #   K<-5
  #   w_train_val<-W_calculus(k = k * n , num_class = K , 
  #                           X = Vect_Unit_Pareto_margins, 
  #                           grid = Grid_points,
  #                           q = q)
  #   Z<-1
  #   for(lam in Vect_lambda){
  #     CV_j<-cross_validation(d , r,  grid = Grid_points,
  #                            lambda = lam , num_col = NULL, start ,
  #                            type = Mod,
  #                            p = p , w = w_train_val , num_class=K)
  #     Optimization_lam_j<-param_estim(d =d, r =r, grid = Grid_points,
  #                                     lambda=lam ,
  #                                     num_col = NULL ,
  #                                     start = start, type = Mod,
  #                                     p = p ,  w = w_total )
  #     L_Lambda[[Z]]<-list("cv"=CV_j,
  #                         "result_theta_A"=Optimization_lam_j)
  #     Z<-Z+1
  #   }
  #   Result_CV<-sapply(L_Lambda,function(x){return(x$cv)})
  #   plot(Vect_lambda,Result_CV,
  #        main="CV score according to lambda value")
  #   Theta_opt<-A_Theta
  # }
  # 
  Approach<-Params_risk_Function[["RF_Approach"]]
  if(Approach=="AD"){
    if(Params_risk_Function[["parametric"]]==TRUE){
      Theta_opt<-Estim_param_RF_homogeneous(Params_risk_Function= Params_risk_Function,
                                            Q_thresh =1-Q_thresh,
                                            Vect_l_function = Vect_l_transf,
                                            d = ncol(Vect_l_transf))
      if(Params_risk_Function[["name_RF"]]=="max"){
        Result_AIC<-as.numeric(as.data.frame(Theta_opt$result_AIC))
        Names_model<-names(Theta_opt$result_AIC)
        Result_AIC<-cbind.data.frame(Names_model,
                                     Result_AIC)
        colnames(Result_AIC)<-c("model","AIC_found")
        Result_AIC$index<-c(1:nrow(Result_AIC))
        GG1<-ggplot(data=Result_AIC,aes(x=index,y=AIC_found,col=model))+
          geom_point()+
          labs(col="MEV model")
        ggsave(filename= paste0(root_for_export,
                                dir_AIC,"_",
                                Approach,".png"),
               plot=GG1,width=8,height=6)
        dev.off()
        Theta_opt<-Theta_opt$alpha
      }
    }
  }
  
  if(Approach=="Gauss"){
    Vect_l_norm<-apply(X = Vect_lunif,
                       FUN = function(x){return(qnorm(x))},
                       MARGIN = 2)
    g_Gauss<-apply(Vect_l_norm,MARGIN = 1,
                   FUN = List_Params_RF[["function"]])
    Mu_Gauss<-quantile(g_Gauss,Q_thresh)
    
    Mu_vector<-colMeans(Vect_l_norm)
    Cov_mat<-cov(Vect_l_norm)
    ### plot MV diag plot
    MVN_diag_plot<-MVN::multivariate_diagnostic_plot(Vect_l_norm)
    png(filename = paste0(root_for_export,"RisKfunctions/",
                          l_name[1],"_",l_name[2],"_Maha_distances.png"))
    print(MVN_diag_plot)
    dev.off()
    ### Summary of results
    Asymp_test<-MVN::mvn(Vect_l_norm, mvn_test="mardia",
                         bootstrap = FALSE)
    print(Asymp_test$multivariate_normality)
    print("### MVN test")
    Mardia_test<-round(MVN::mardia(Vect_l_norm)$p.value,2)
    print(Mardia_test)
    Royston_test<-round(MVN::royston(Vect_l_norm)$p.value,2)
    Energy_test<-round(MVN::energy(Vect_l_norm)$p.value,2)
    HZ_test<-round(MVN::hz(Vect_l_norm)$p.value,2)
    Name_models<-c("mardiaSkew","mardiaKurtosis","Royston",
                   "Energy test","HZ")
    Result_models<-c(Mardia_test, Royston_test,
                     Energy_test,HZ_test )
    Table_Gauss_hyp<-cbind.data.frame(Name_models,Result_models)
    colnames(Table_Gauss_hyp)<-c("name_test",
                                 "p_value")
    write.csv(file=paste0(root_for_export,"RiskFunctions/",
                          l_name[1],"_",l_name[2],"_MVN_test_results.csv"),
              x = Table_Gauss_hyp)
  }
  if(Approach=="HTawn"){
    
    # Vect_l_marg<-apply(X = Vect_lunif,
    #                    FUN = function(x){return(qnorm(x))},
    #                    MARGIN = 2)
    # g_marg<-apply(Vect_l_marg,MARGIN = 1,
    #               FUN = List_Params_RF[["function"]])
    # Mu_marg<-quantile(g_marg,Q_thresh)
    Vect_log<-log(Vect_l_function+1)
    g_marg<-apply(Vect_log,MARGIN = 1,
                  FUN = List_Params_RF[["function"]])
    Mu_marg<-quantile(g_marg,
                      Q_thresh)
    IREF<-1
    ### probability of being in each case
    #Vect_l_marg
    Q_from_mixt<-apply(Vect_log,
                       MARGIN = 2,FUN = function(x){
                         return(mean(as.numeric(x>Mu_marg)))})
    DQU_each_scenario<-1-Q_from_mixt
    ### Choose the du used during modelling step.
    # Use the one used for future simulations
    # DQU_modeling_Htawn<-Q_from_mixt
    d<-length(l_variables)
    DQU_modeling_Htawn<-rep(NA,d)
    nb_threshs<-50
    Qtile_candidates<-seq.int(0.65,0.85,
                              length.out=nb_threshs)
    MQU<-0.70
    DQU_modeling_Htawn<-Params_risk_Function[["HTAWN_params"]]
    print(DQU_modeling_Htawn)
    Graphics_diags<-Analysis_diag_HTawn_evol_DQU(vect_dqu = Qtile_candidates,
                                                 chosen_dqu=DQU_modeling_Htawn,
                                                 MQU =MQU,
                                                 Vect_obs = Vect_log)
    GG_theta<-Graphics_diags[["theta"]]
    GG_Indep<-Graphics_diags[["indep"]]
    ggsave(filename = paste0(root_for_export,"RiskFunctions/",
                             l_name[1],"_",l_name[2],
                             "_evol_HTparams.png"),
           plot = GG_theta,width = 8,height = 6)
    ggsave(filename = paste0(root_for_export,"RiskFunctions/",
                             l_name[1],"_",l_name[2],
                             "_result_HTindep.png"),
           plot =  GG_Indep,width = 8,height = 6)
    ### HTawn modelling
    # vect_l_marg
    Model_Htawn_all<-texmex::mexAll(Vect_log,
                                    mqu = MQU,
                                    dqu = DQU_modeling_Htawn)
    Theta_opt<-Model_Htawn_all
    list_GG_cases<-list()
    for(Z in c(1:d)){
      HTawn_z<-Theta_opt[[Z]]
      ### Graphics to visualize what appears
      Data_tfed<-as.data.frame(HTawn_z$margins$transformed[,l_variables])
      colnames(Data_tfed)<-sapply(c(1:d),function(x){
        return(paste0("V",x))
      })
      Th_used_j<-as.numeric(HTawn_z$dependence$dth)
      Data_tfed$inds_exts_Z<-as.character(
        as.numeric(Data_tfed[,Z]>Th_used_j))
      GG_z_model<-ggplot(data = Data_tfed,aes(x=V1,y=V2,
                                              col=inds_exts_Z))+
        geom_point()
      if(Z==1){
        GG_z_model<-GG_z_model+
          geom_vline(aes(xintercept=Th_used_j))
      }else{
        GG_z_model<-GG_z_model+
          geom_hline(aes(yintercept=Th_used_j))
      }
      GG_z_model<-GG_z_model+ 
        labs(col="Legend")+
        xlab(paste0("Norm of T(",l_variables[1],")"))+
        ylab(paste0("Norm of T(",l_variables[2],")"))
      list_GG_cases[[Z]]<-GG_z_model
    }
    GG_hull<-Analysis_diag_HTawn_Conv_Hull(
      Theta_opt = Theta_opt)
    ggsave(plot = GG_hull,filename = paste0(root_for_export,"RiskFunctions/",
                                            l_name[1],"_",l_name[2],
                                            "_Convex_huls_Htawn.png"),
           width=8,height=8)
    
    ### RiskFunctions
    Object_cases_HTAWN<-list_GG_cases[[1]]+list_GG_cases[[2]]+
      plot_layout(axis_titles = "collect") 
    
    ggsave(plot = Object_cases_HTAWN,
           filename = paste0(root_for_export,"RiskFunctions/",
                             l_name[1],"_",l_name[2],
                             "_pts_input_HT.png"),
           width=8,height=6)
    ### Simulation step
    prop_cl<-Q_from_mixt
    # AI case
    # directly Q_from_mixt
    # AD case) other j
    W_scenarios<-rep(Qfound,d)
    W_scenarios[IREF]<-Q_from_mixt[IREF]
    W_scenarios[-IREF]<-W_scenarios[-IREF]-(Q_from_mixt[IREF])
    probs_cl<-W_scenarios/sum(W_scenarios)
  }
  if(Approach=="Laplace"){
    Vect_Laplace<-apply(X = Vect_lunif,
                        FUN = function(x){
                          return(L1pack::qlaplace(x))},
                        MARGIN = 2)
    g_Laplace<-apply(Vect_Laplace,MARGIN = 1,
                     FUN = List_Params_RF[["function"]])
    Mu_Laplace<-as.numeric(quantile(g_Laplace,Q_thresh))
    Model_ALD<-L1pack::LaplaceFit(x = Vect_Laplace)
  }
  if(Approach=="HGD"){
    d<-ncol(Vect_l_norm)
    model_univariate<-ghyp::ghyp()
    Norm_std<-apply(Vect_lunif,MARGIN=2,FUN = Convert_marg,
                    model_univ=model_univariate)
    Model_HGD<-ghyp::fit.VGmv(data = Norm_std,
                              gamma=rep(0,d),
                              lambda = 1,opt.pars = c(lambda=TRUE,
                                                      mu=TRUE,gamma=FALSE,
                                                      sigma=TRUE),
                              symmetric = TRUE,
                              sigma=cov(Norm_std),
                              mu=rep(0,d),
                              standardize=TRUE)
    # Model_HGD<-Fit_HGD_from_Laplace(Lap_vectors = Vect_Laplace,
    #           Model_Lap_inits = Model_ALD)
  }
  MATRICE_SCORES<-c()
  l_fonction<-list()
  l_Mat_moyenne<-list()
  l_sigma<-list()
  LISTE_SHAPE_OBS<-list()
  d<-ncol(Vect_l_function)
  
  L<-length(l_variables)-1
  for(NameVAR in l_variables[1:L]){
    Name_for_export<-paste0(root_for_export,NameVAR,"_")
  }
  Name_for_export<-paste0(Name_for_export,"_",l_variables[length(l_variables)])
  if(One_PCA_base){
    print("One PCA basis")
    LIST_Mod<-Approach_Angle_One_PCA(LIST_all = Liste_all,
                                     Name_for_export =root_for_export,
                                     NbScores_Omega = NbScores_Omega,
                                     l_variables = l_variables,
                                     Indices_exts = Indices_exts,
                                     d = d,f_transf = f_transf)
  }else{
    print("Several PCA basis")
    LIST_Mod<-Approach_Angle_Mult_PCA(Indices_exts = Indices_exts,
                                      LIST_all = Liste_all,
                                      root_export = root_for_export,
                                      l_variables = l_variables,
                                      list_nb_scores=list_nb_scores,
                                      f_transf = f_transf)
    NbScores_Omega<-LIST_Mod[["Nb_scores_omega"]]
  }
  LISTE_Frechet_OBS<-LIST_Mod[["L_Frechet"]]
  Scores<-LIST_Mod[["Scores"]]
  Length_T<-LIST_Mod[["Length_TS"]]
  EIGEN_functions<-LIST_Mod[["Eigen_functions"]]
  Unif_coord<-VineCopula::pobs(Scores)
  M_Theta<-Params_risk_Function[["Method_Angle"]]
  if(M_Theta=="VineCop"){
    Matrix_C<-function_Structure_Matrice(NB_dim = NbScores_Omega)
    Modele_coords<-VineCopula::RVineCopSelect(Unif_coord,
                                              Matrix = Matrix_C)
  }
  if(M_Theta=="BetaCop"){
    Modele_coords<-copula::empCopula(X = Unif_coord,
                                     smoothing = "beta")
  }
  if(M_Theta=="GaussMixture"){
    Mod_Gauss<-mclust::densityMclust(Scores,
                                     plot=FALSE)
  }
  Bool_filled<-FALSE
  N_sim<-M
  Coords_simul<-matrix(NA,nrow = M,ncol = NbScores_Omega)
  while(Bool_filled==FALSE){
    if(M_Theta!="GaussMixture"){
      if(M_Theta=="VineCop"){
        # Simulations of coordinates  ------------------------------------------------------------
        Simulations_coord<-VineCopula::RVineSim(N=N_sim,
                                                RVM = Modele_coords)
      }
      if(M_Theta=="BetaCop"){
        Simulations_coord<-copula::rCopula(N_sim, 
                                           copula=Modele_coords)
      }
      
      colnames(Simulations_coord)<-1:ncol(Simulations_coord)
      Coords_ech_orig<-matrix(NA,nrow=N_sim,
                              ncol=NbScores_Omega)
      for(j in c(1:ncol(Simulations_coord))){
        sortie<-quantile(Scores[,j],Simulations_coord[,j])
        Coords_ech_orig[,j]<-as.numeric(sortie)
      }
    }
    if(M_Theta=="GaussMixture"){
      Coords_ech_orig<-t(replicate(n = N_sim,
                                   expr=Simul_from_MixtureGauss(Mod_Gauss))[1,,])
    }
    if(One_PCA_base){
      Shape_Omega_simul<-Simul_Omega_One_PCA_base(list_Mod_One_PCA = LIST_Mod,
                                                  Simul_coords =Coords_ech_orig 
                                                  ,NbScores_Omega =NbScores_Omega,
                                                  f_transf_inv = f_transf_inv)
    }else{
      Shape_Omega_simul<-Simul_Omega_Mult_PCA_base(list_Mod_Mult_PCA =LIST_Mod,
                                                   M = N_sim,d = d,list_nb_scores = list_nb_scores,
                                                   Simul_coords =Coords_ech_orig,
                                                   f_transf_inv = f_transf_inv)
    }
    
    ##### Compar simul of theta with the reality
    scale_frechet_d<-Params_risk_Function[["scale_frechet_d"]]
    if(Approach=="AD"){
      if(Params_risk_Function[["parametric"]]==FALSE){
        UNIF<-exp(-Vect_l_transf^(-1))
        EXP<--log(1-UNIF)
        Sim_l<-Simul_RF_max_NP(Params_risk_Function = Params_risk_Function,
                               M=N_sim,QTH = Q_thresh,
                               Vect_l_function = EXP)
        Lg<-apply(X =EXP,MARGIN = 1,FUN = max)
        Threshold_exp<-quantile(Lg,Q_thresh)
        Sim_l_tf<-pexp(Sim_l+Threshold_exp)
      }else{ 
        Sim_l<-Simulation_gParetoP(Params_risk_Function=Params_risk_Function,
                                   d=ncol(Vect_l_function),
                                   M=N_sim,theta_opt=Theta_opt)
        ### From Legrand
        ### Warning) definition by continuity if null values are given
        Sim_l<-Sim_l*Seuil_Margins_Pareto
        Sim_l_tf<-apply(Sim_l,MARGIN=2,FUN = evd::pgpd,loc=0,
                        scale=1, shape=1)
      }
    }
    if(Approach=="Gauss"){
      ### To do: correct Sample_cond_g
      Sim_l<-t(replicate(n =N_sim,
                         expr = Sample_cond_g_Rjection_Sampling(Mu_vector = Mu_vector,
                                                                Cov_mat = Cov_mat,
                                                                g=List_Params_RF[["function"]],
                                                                Th_g=Mu_Gauss)))
      ### Conversion to uniform margins.
      print("here")
      Sim_l_tf<-apply(Sim_l,MARGIN=2,
                      FUN = pnorm)
      ### Mixture model
      
    }
    if(Approach=="Laplace"){
      Sim_l<-t(replicate(n =N_sim,
                         expr = Sample_cond_Laplacian_Rjection_Sampling(
                           Mu_ALD = Model_ALD$center,
                           Scatter_ALD=Model_ALD$Scatter,
                           g=List_Params_RF[["function"]],
                           Th_g=Mu_Laplace)))
      Sim_l_tf<-apply(Sim_l,MARGIN=2,
                      FUN = L1pack::plaplace)
    }
    if(Approach=="HGD"){
      Sim_l<-t(replicate(n =N_sim,
                         expr = Sample_cond_HGD_Rjection_Sampling(
                           Model_HGD = Model_HGD,
                           g=List_Params_RF[["function"]],
                           Th_g=Mu_Gauss)))
      Sim_l_tf<-apply(Sim_l,MARGIN=2,
                      FUN =pnorm)
    }
    if(Approach=="HTawn"){
      Sim_l<-Simul_from_Htawn(model_mex_all = Model_Htawn_all,
                              n_sim = N_sim,d = d,
                              Name_vars = l_variables,
                              prop_class = prop_cl,ind_ref=IREF,
                              thresh_censor=as.numeric(Mu_marg),
                              Q_sim = DQU_each_scenario)
      # Sim_l_tf<-apply(Sim_l,MARGIN=2,
      #                 FUN = pnorm)
      Sim_l<-exp(Sim_l)-1
    }
    
    if(Approach!="HTawn"){
      ### We cannot use the HTawn modelling of the radial
      ### components.
      Sim_l<-Sim_l_tf
      for(j in c(1:length(l_variables))){
        var_tf_sim<-Sim_l_tf[,j]
        nameV<-l_variables[j]
        Param_tf_lj<-L_Transf_vect_l[[nameV]]
        # WEIGHT<-Param_tf_lj[["p_u"]]
        # R_j<-sapply(var_tf_sim,Fast_reconv,p_Ui=WEIGHT,
        #             scale_fonc=Param_tf_lj[["scale"]],
        #             gamma=Param_tf_lj[["gamma"]],
        #             kernel_DENS=Param_tf_lj[["kernel_dens"]],
        #             u_f=Param_tf_lj[["threshold"]])
        #ELT_pos<-which(var_tf_sim>0)
        R_j<-var_tf_sim
        R_j<-sapply(var_tf_sim,mev::qextgp,
                    kappa =Param_tf_lj[["kappa"]],
                    sigma = Param_tf_lj[["sigma"]],
                    xi = Param_tf_lj[["xi"]])  
        Sim_l[,j]<-exp(R_j)-1
      }
    }
    
    ### Extreme individuals in the original scale Pareto(shape param)
    Excedents_lprime_orig<-Vect_l_function[Indices_exts,]
    plot(Excedents_lprime_orig,log="xy")
    points(Sim_l,col="blue",
           cex=0.75)
    print("Kendall correlation")
    Kendall_sim<-cor(Sim_l[,1],Sim_l[,2],
                     method = "kendall")
    print(Kendall_sim)
    Kendall_obs<-cor(Excedents_lprime_orig[,1],
                     Excedents_lprime_orig[,2],
                     method = "kendall")
    print(Kendall_obs)
    ### Chimeas for sim L (vs) obs L
    png(filename= paste0(root_for_export,"RisKfunctions/chi_meas_sim_obs",
                         Params_risk_Function[["general_option"]],".png"),
        width=750,250)
    par(mfrow=c(1,2))
    POT::chimeas(Excedents_lprime_orig,which=1)
    POT::chimeas(Sim_l,which=2)
    par(mfrow=c(1,1))
    dev.off()
    ####
    Sim_l<-as.data.frame(Sim_l)
    colnames(Sim_l)<-l_variables
    Df_combs_sim_exts<-rbind.data.frame(Excedents_lprime_orig,Sim_l)
    colnames(Df_combs_sim_exts)<-c("V1","V2")
    Df_combs_sim_exts$Legend<-c(rep("data",
                                    nrow(Excedents_lprime_orig)),
                                rep("simulations",
                                    nrow(Sim_l)))
    Z_gg<-1
    Z_gg2<-2
    XLAB<-expression(R[M~","~Z_gg]^g)
    XLAB<-as.expression(do.call('substitute', list( XLAB[[1]], 
                                                    list(Z_gg=Z_gg))))
    YLAB<-expression(R[M~","~Z_gg2]^g)
    YLAB<-as.expression(do.call('substitute', 
                                list( YLAB[[1]], list(Z_gg2=Z_gg2))))
    GG_RG_sim_vs_obs<-ggplot(data=Df_combs_sim_exts,
                             aes(x=V1,y=V2))+
      geom_point(aes(shape=Legend,col=Legend,
                     size=Legend))+
      geom_xsidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                        alpha = 0.5)+
      geom_ysidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                        alpha = 0.5)+
      scale_y_continuous(transform = "log10")+
      scale_x_continuous(transform = "log10")+
      # xlab(paste0("l(T(",l_variables[1],"))"))+
      # ylab(paste0("l(T(",l_variables[2],"))"))+
      xlab(XLAB)+
      ylab(YLAB)+
      scale_shape_manual(values = c("simulations"=17,
                                    "data"=19))+
      scale_size_manual(values=c("simulations"=0.75,
                                 "data"=1.5))+
      scale_color_manual(values=cols_ggplot)+
      scale_fill_manual(values=cols_ggplot)+
      guides(fill="none")+
      theme(axis.title=element_text(size=15),
            legend.text=element_text(size=10))
    ggsave(filename= paste0(root_for_export,"RisKfunctions/Rg_sim_vs_obs_",
                            l_variables[1],"_",l_variables[2],"_",
                            Params_risk_Function[["general_option"]],".png"),
           plot=GG_RG_sim_vs_obs,
           width=8,height=6)
    LISTE_shapes<-list()
    vect_found<-c()
    LISTE_candidats<-list()
    DF_simul_lprime<-c()
    Z<-1
    for(k in c(1:length(l_variables))){
      nom_s<-l_variables[k]
      L_prime_simul<-Sim_l[,k]
      #Deconcatenate the Omega----------------
      #####
      Beg<-1+Length_T*(Z-1)
      END<-Length_T*Z
      Shape_forcing<-Shape_Omega_simul[,Beg:END]
      L2_shape_name<-apply(Shape_forcing,MARGIN = 1,FUN = calcul_norme_L2)
      Shape_forcing_std<-t(t(Shape_forcing)%*%diag(L2_shape_name^(-1)))
      LISTE_shapes[[nom_s]]<-Shape_forcing_std
      Z_varj<-t(t(Shape_forcing_std)%*%(diag(L_prime_simul)))
      LISTE_candidats[[nom_s]]<-Z_varj
      DF_simul_lprime<-c(DF_simul_lprime,
                         apply(X = Z_varj,MARGIN = 1,
                               FUN = calcul_norme_L2))
      ### 0 as the threshold if Frechet
      ### 1 as the threshold if Pareto
      ### If Frechet margins originally
      #Unif<-exp(-Z_varj[indicatrice_pos,]^(-1))
      ### If Pareto margins originally
      # if(opt_Frech){
      #   indicatrice_pos<-which(apply(X=Z_varj,
      #                                FUN=fonction_trajectoire_positive,
      #                                MARGIN = 1)==TRUE)
      #   Unif<-exp(-(Z_varj[indicatrice_pos,])^(-1))
      # }else{
      #   indicatrice_pos<-which(apply(X=Z_varj-1,
      #                                FUN=fonction_trajectoire_positive,
      #                                MARGIN = 1)==TRUE)
      #   Unif<-1-Z_varj[indicatrice_pos,]^(-1)
      # }
      indicatrice_pos<-which(apply(X=Z_varj,
                                   FUN=fonction_trajectoire_positive,
                                   MARGIN = 1)==TRUE)
      Unif<-apply(Z_varj[indicatrice_pos,],
                  MARGIN=2,FUN = evd::pgpd,loc=0,
                  scale=Shape_evd, shape=Shape_evd)
      PTO<-(1-Unif)^(-1)
      seuil_marg<-1+10^(-5)
      indicatrice_pos2<-which(apply(X=(PTO-seuil_marg),
                                    FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
      indicatrice_pos<-indicatrice_pos[indicatrice_pos2]
      vect_found<-c(vect_found,indicatrice_pos)
      Z<-Z+1
    }
    ref<-Array_simul[1,,1]
    Nb_filled<-length(which(!is.na(ref)==TRUE))
    Inds_final_chosen<-unique(vect_found[duplicated(vect_found)])
    Nb_final_chosen<-length(Inds_final_chosen)
    Beg_filled<-Nb_filled+1
    Gap<-M-(Nb_filled+Nb_final_chosen)
    End_filled<-min(Nb_filled+Nb_final_chosen,M)
    Lag<-End_filled-Beg_filled+1
    for(k in c(1:length(l_variables))){
      name_variable<-l_variables[k]
      Z_vark<-LISTE_candidats[[name_variable]][Inds_final_chosen,]
      Sub_coords<-Coords_ech_orig[Inds_final_chosen,]
      Sub_coords<-Sub_coords[c(1:Lag),]
      # UNIF_k<-exp(-Z_vark^(-1))
      # if(opt_Frech){
      #   
      # }else{
      #   UNIF_k<-1-Z_vark^(-1)
      # }
      UNIF_k<-apply(Z_vark,MARGIN=2,
                    FUN = evd::pgpd,loc=0,
                    scale=Shape_evd, shape=Shape_evd)
      UNIF_k<-UNIF_k[c(1:Lag),]
      ### Unif scale
      Array_simul[k,c(Beg_filled:End_filled),]<-UNIF_k
      Coords_simul[c(Beg_filled:End_filled),]<-Sub_coords
    }
    ### See if the object is completed. 
    new_ref<-Array_simul[1,,1]
    Nb_filled_update<-length(which(!is.na(new_ref)==TRUE))
    ### if True, stop the loop
    if(Nb_filled_update==M){
      Bool_filled<-TRUE
    }
  }
  
  # Conversion --------------------------------------------------------------
  # ------------------------------------------------------------------------
  LISTE_Frechet_SIMUL<-list()
  LISTE_simul<-list()
  LISTE_obs_exts<-list()
  if(result_transformation[["type_transfo"]]=="Mixture_emp_GPD"){
    for(k in c(1:length(l_variables))){
      name_variable_for_conv<-l_variables[k]
      Simul_whole_kunif<-Array_simul[k,,]
      ### Unif-->Pareto to compare in the scale of obs.
      # If GPD(1,0,1) in conversion
      Conv_for_analyse<-apply(Simul_whole_kunif,MARGIN=2,
                              FUN = evd::qgpd,loc=0,
                              scale=Shape_evd, shape=Shape_evd)
      # if Frechet in conversion. 
      #Conv_for_analyse<-(1-Simul_whole_kunif)**(-1)
      LISTE_Frechet_SIMUL[[name_variable_for_conv]]<-Conv_for_analyse
      INDICES_ACP<-1:nrow(Z_varj)
      ### Unif--> Pareto for conversion
      NO_ACP<-lapply(INDICES_ACP,FUN = fnct_select_colonne,
                     df=(1-Simul_whole_kunif^(-1)))
      
      # Reconversion in the good scale--------------------------------------
      K<-Liste_all[[name_variable_for_conv]]$K
      LEVT<-Liste_all[[name_variable_for_conv]]$LEVT
      Variables_reconversion_ACP<-lapply(NO_ACP,function_reconversion_Pareto,
                                         K=K,
                                         list_evt=LEVT)
      Variables_reconversion_ACP<-do.call(rbind.data.frame,
                                          Variables_reconversion_ACP)
      colnames(Variables_reconversion_ACP)<-1:ncol(Variables_reconversion_ACP)
      LISTE_simul[[name_variable_for_conv]]<-Variables_reconversion_ACP
      
      # Observation extreme -----------------------------------------------------
      obs_ext<-result_transformation$orig[[name_variable_for_conv]][Indices_exts,]
      colnames(obs_ext)<-c(1:ncol(obs_ext))
      LISTE_obs_exts[[name_variable_for_conv]]<-obs_ext
    }
  }else{
    for(k in c(1:length(l_variables))){
      name_variable_for_conv<-l_variables[k]
      Simul_whole_kunif<-Array_simul[k,,]
      #Conv_for_analyse<-(-1)*(log(Simul_whole_kunif))^(-1)
      ### Unif-->Pareto to compare in the scale of obs.
      Conv_for_analyse<-apply(Simul_whole_kunif,MARGIN=2,
                              FUN = evd::qgpd,loc=0,
                              scale=Shape_evd, 
                              shape=Shape_evd)
      #Conv_for_analyse<-(1-Simul_whole_kunif)**(-1)
      INDICES_ACP<-1:nrow(Simul_whole_kunif)
      ### Unif-->Frechet to compare in Frechet scale of obs.
      LISTE_Frechet_SIMUL[[name_variable_for_conv]]<-Conv_for_analyse
      NO_ACP<-lapply(INDICES_ACP,FUN = fnct_select_colonne,
                     df=Simul_whole_kunif)
      
      # Reconversion in the good scale--------------------------------------
      Theta_EXTGPD_k<-Liste_all[[name_variable_for_conv]]$params_transfo
      Variables_reconversion_ACP<-lapply(NO_ACP,
                                         function_reconversion_Unif_EXTGPD,
                                         Theta_k=Theta_EXTGPD_k)
      Variables_reconversion_ACP<-do.call(rbind.data.frame,
                                          Variables_reconversion_ACP)
      colnames(Variables_reconversion_ACP)<-1:ncol(Variables_reconversion_ACP)
      LISTE_simul[[name_variable_for_conv]]<-Variables_reconversion_ACP
      
      # Observation extreme -----------------------------------------------------
      obs_ext<-result_transformation$orig[[name_variable_for_conv]][Indices_exts,]
      colnames(obs_ext)<-c(1:ncol(obs_ext))
      LISTE_obs_exts[[name_variable_for_conv]]<-obs_ext
    }
    
  }
  if(M_Theta=="GaussMixture"){
    Object_modelisation<-Mod_Gauss
  }
  else{
    Object_modelisation<-Modele_coords
  }
  return(list("Risk_MV"=l_prime[Indices_exts],
              "Param_found"=Theta_opt,"simul"=LISTE_simul,
              "obs_exts"=LISTE_obs_exts,"Indices_exts"=Indices_exts,
              "coords_simul"=Coords_simul, "coords_data"=Scores,
              "Frechet_normal"=list("obs"=LISTE_Frechet_OBS,
                                    "simul"=LISTE_Frechet_SIMUL),
              "EIGEN_functions"= EIGEN_functions,
              "object_model_theta"=Object_modelisation,
              "threshold"=Seuil_lprime,
              "Vect_RiskF_transf"=Vect_l_transf,
              "model_EGPD"=L_Transf_vect_l))
}

#' function_reconversion_Unif_EXTGPD
#'
#' @param Simul_Unif: vector. Simulated vector with uniform margins.
#' @param Theta_k: list. List of EGPD estimators found at each time t.
#'
#' @return Simulated vector with EGPD margins.
#' @export
#'
#' @examples
function_reconversion_Unif_EXTGPD<-function(Simul_Unif,Theta_k){
  
  L_dim<-length(Simul_Unif)
  individu_ext<-sapply(c(1:L_dim),function_reconv_each_dim_EGPD,Simul_Unif=Simul_Unif,
                       Theta_k=Theta_k)
  return(individu_ext)
}
#' function_reconv_each_dim_EGPD
#'
#' @param Simul_Unif: vector. Simulated vector with uniform margins.
#' @param Theta_k: list. List of EGPD estimators found at each time t.
#' @param index_dim: int. Time index. 
#'
#' @return Value at time t of the vector at EGPD scale. 
#' @export
#'
#' @examples
function_reconv_each_dim_EGPD<-function(Simul_Unif,Theta_k,index_dim){
  
  value_unif_x<-Simul_Unif[index_dim]
  Theta_kt<-Theta_k[[index_dim]]
  Nu_t<-Theta_kt[["nu"]]
  Gamma_t<-Theta_kt[["mu"]]
  Sigma_t<-Theta_kt[["sigma"]]
  #Apply inverse of G
  value_G_x<-value_unif_x^(1/Nu_t)
  #Apply inverse of GPD cd
  value_EGPD_t<-((1-value_G_x)^(-Gamma_t)-1)*(Sigma_t/Gamma_t)
  return(value_EGPD_t)
}


Function_one_couple_l<-function(Params_Biv,z,Th,one_u){
  
  realisation<-Function_mixture(Params_Biv =Params_Biv ,u = one_u,
                                Th =Th,z = z)
  realisation_y<-Simul_cond_one(input_x = realisation,
                                Params_biv = Params_Biv,
                                Th = Th,
                                z = z)
  return(c(realisation,realisation_y))
}

Min_function_conv<-function(Shape_d){
  pas_x<-1/ncol(Shape_d)
  vecteur_temps<-c(1:ncol(Shape_d))/ncol(Shape_d)
  #Fonction propre. 
  fonction_propre_j<-function(vecteur_temps,j){
    fnct_par_temps<-function(j,t){
      return(sin(2*pi*t*j))
    }
    vecteur_r<-sapply(vecteur_temps,fnct_par_temps,j=j)
    return(vecteur_r)
  }
  LISTE_convergence<-c()
  for(l in c(1:8)){
    fonction_obtenue<-fonction_propre_j(vecteur_temps =vecteur_temps,j=l )
    
    # approx de Rieman --------------------------------------------------------
    coordonnees<-(Shape_d%*%fonction_obtenue)*(pas_x)
    LISTE_convergence<-c(LISTE_convergence,mean(abs(coordonnees)))
  }
  
  return(LISTE_convergence)
}

DIndex_Poisson_Thresholds<-function(Vect_l_function,U,V,Decimal_date,Nb_years){
  NPobs<-nrow(Vect_l_function)/Nb_years
  Result<-1-as.numeric((Vect_l_function[,1]/NPobs<U)&(Vect_l_function[,2]/NPobs<V))
  Lambda<-sum(Result)/Nb_years
  #from diplot
  date <- floor(Decimal_date)
  tim.rec <- range(date)
  nb.occ <- NULL
  for (year in tim.rec[1]:tim.rec[2]) nb.occ <- c(nb.occ, 
                                                  sum(Result & (date == year)))
  Variance_<-var(nb.occ)
  DIndex<-Lambda/Variance_
  return(DIndex)
}
Analyse_extreme_proj_gfunction<-function(base_RV,L,M1,M2,function_g){
  
  fnct_k<-function(Obs,k,function_g){
    LISTE_p<-function_analyse_convergence_gPto(Obs=Obs,K= k,
                                               function_g=function_g)
    return(LISTE_p)
  }
  vecteur_k<-seq(50,L,by=1)
  # Travail sur la convergence vers un processus l-Pareto. ----------------------------------------------
  
  vecteur_CONVERG<-sapply(X = vecteur_k,FUN = fnct_k,Obs=base_RV,
                          function_g=function_g)
  MATRICE_moy_coord<-t(cbind.data.frame(vecteur_CONVERG))
  par(mfrow=c(2,2))
  # TITLE_proj<-expression("First moment of "~Theta[M]~"'s projection for "~nom_graph)
  # TITLE_proj<-do.call("substitute", list(TITLE_proj[[1]], list(nom_graph = nom_graph)))                                                                            
  for(j in c(1:ncol(MATRICE_moy_coord))){
    express_h<-expression("Moment with "~h[j])
    express_h<-do.call("substitute", list(express_h[[1]], list(j = j)))
    plot(MATRICE_moy_coord[,j],type="l",ylab=express_h,
         xlab="Exceedances",cex.lab=1.2)
    abline(v=M1,col="red")
    abline(v=M2,col="red")
    # if(j==4){
    #   mtext(outer=TRUE,text =  TITLE_proj,
    #         line = -2)
    # }
  }
  # mtext(outer=TRUE,text = TITLE_proj,
  #       line = -2)
  par(mfrow=c(1,1))
}

function_analyse_convergence_gPto<-function(Obs,K,function_g){
  
  g_data<-apply(Obs,FUN =function_g,MARGIN = 1)
  pas_x<-1/ncol(Obs)
  indice_k<-order(g_data,decreasing = TRUE)[K]
  Seuil_L2<-g_data[indice_k]
  Indices<-which(g_data>Seuil_L2)
  Conservees<-Obs[Indices,]
  g_data_cons<-g_data[Indices]
  FORME_d<-t(t(Conservees)%*%diag(g_data_cons^(-1)))
  vecteur_temps<-c(1:ncol(Obs))/ncol(Obs)
  #Fonction propre. 
  fonction_propre_j<-function(vecteur_temps,j){
    fnct_par_temps<-function(j,t){
      return(sin(2*pi*t*j))
    }
    vecteur_r<-sapply(vecteur_temps,fnct_par_temps,j=j)
    return(vecteur_r)
  }
  LISTE_convergence<-c()
  for(l in c(1:8)){
    fonction_obtenue<-fonction_propre_j(vecteur_temps =vecteur_temps,j=l )
    
    # approx de Rieman --------------------------------------------------------
    coordonnees<-(FORME_d%*%fonction_obtenue)*(pas_x)
    LISTE_convergence<-c(LISTE_convergence,mean(abs(coordonnees)))
  }
  
  return(LISTE_convergence)
}

### base comes from tea::mindist, idea is to better 
### see the variations of the distance 
### to prevent take a not so interesting threshold.
mindist_update<-function (data, ts = 0.15, method = "mad") 
{
  xstat = sort(data, decreasing = TRUE)
  n = length(data)
  T = floor(n * ts)
  i = 1:(n - 1)
  h = (cumsum(log(xstat[i]))/i) - log(xstat[i + 1])
  xstat = sort(data)
  A = matrix(ncol = T - 1, nrow = T - 1)
  for (k in 1:(T - 1)) {
    for (j in 1:(T - 1)) {
      A[k, j] = abs((((k/j) * xstat[n - k + 1]^(1/h[k]))^h[k]) - 
                      xstat[n - j])
    }
  }
  if (method == "mad") {
    M = rowMeans(A)
  }
  if (method == "ks") {
    rowMax <- function(rowData) {
      apply(rowData, MARGIN = c(1), max)
    }
    M = rowMax(A)
  }
  method<-rep(method,length(M))
  df<-cbind.data.frame(1:(T-1),M,method)
  colnames(df)<-c("Nb_k","value_metric","metric")
  return(df)
}
#' Couple_s_t_chi_measure
#'
#' @param t: int. First time of the pair.  
#' @param s: int. Second time of the pair. 
#' @param matrix_: matrix. Observations matrix.
#' @param tail_quantile: float. High threshold in the uniform scale. 
#'
#' @return Extremal correlation coefficient Chi for the pair at (t,s) 
#' @export
#'
#' @examples
Couple_s_t_chi_measure<-function(t,s,matrix_,tail_quantile){
  
  Sub_matrix<-matrix_[,c(t,s)]
  Min_st<-apply(X = Sub_matrix,MARGIN = 1,FUN = min)
  prob_intersection<-mean(as.numeric(Min_st>tail_quantile))
  return(prob_intersection/(1-tail_quantile))
}
#'Matrix_chi_measure
#'
#' @param Data_unif: matrix. Observation with uniform margins. 
#' @param tail_quantile: float. High threshold in the uniform scale. 
#'
#' @return Matrix of extremal correlation coefficient Chi per pair (t,s) 
#' @export
#'
#' @examples
Matrix_chi_measure<-function(Data_unif,tail_quantile){
  
  matrix_chi<-matrix(0,nrow = ncol(Data_unif),
                     ncol=ncol(Data_unif))
  for(i in c(1:ncol(matrix_chi))){
    sup_ij<-c(i:d)
    Values_chi<-sapply(X = sup_ij,FUN = Couple_s_t_chi_measure,
                       t=i,matrix_=Data_unif,tail_quantile=tail_quantile)
    matrix_chi[i,sup_ij]<-Values_chi
  }
  Comb<-t(matrix_chi)+matrix_chi
  diag(Comb)<-diag(Comb)/2
  return(Comb)
}

