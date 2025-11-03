
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
MarTransfo_TS_exts_MV<-function(type_donnees,coeurs,lien_racine,liste_noms,file_dates,
                                p_U,n.dens,opt_Frech,type_entree){
  l_ALL<-list()
  l_Orig<-list()
  j<-1
  l_AD<-list()
  for(nom_variable in liste_noms){
    lien_donnees<-paste0(lien_racine,nom_variable,"_residus.csv")
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
    l_Orig[[nom_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,fonction_seuil = p_U[[nom_variable]],n.dens = n.dens,
                              nom=nom_variable,type_entree=type_entree,
                              dates_prises=dates_import,j_show = 19)
    l_AD[[nom_variable]]<-AN_GPD
    if(type_entree=="donnees_brutes"){
      return(TRUE)
    }
    Vecteur_DIM<-c(1:ncol(Donnes))
    Resultat_P<-as.data.frame(t(sapply(Vecteur_DIM,f_marginales_all_Pareto,data_to_tf=Donnes,
                                       p_u=p_U[[nom_variable]],
                                       n.dens=n.dens)))
    K<-Resultat_P$kernel_dens
    PARETO<-Resultat_P$obs
    L_EVT<-list("gamma"=Resultat_P$gamma,
                "scale"=Resultat_P$scale,
                "threshold"=Resultat_P$threshold,
                "p_u"=Resultat_P$p_u)
    Vtransf<-do.call(cbind.data.frame,PARETO)
    if(opt_Frech==TRUE){
      Vunif<-(1-1/Vtransf)
      Vtransf<--1/log(Vunif)
    }
    l_ALL[[nom_variable]]<-list("transf"=Vtransf,"K"=K,
                                "LEVT"=L_EVT)
    if(nom_variable==liste_noms[1]){
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
  FIT_pos<-mev::gp.fit(xdat = Col_pos_t,threshold =Th)$est
  Shape<-FIT_pos[2]
  # ## Use gpd property to get the scale at 0.
  Scale<-FIT_pos[1]-Th*Shape
  INIT<-c(Shape,Scale)
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
  muFit <- fitted(Fitting_time_t,"mu")[1]
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
              "Init"=INIT))
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
  for(nom_variable in liste_noms){
    lien_donnees<-paste0(lien_racine,nom_variable,"_residuals.csv")
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
  
    l_Orig_all[[nom_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    dates_import<-dates_import[-1]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    print(dim(Donnes))
    AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,fonction_seuil = p_U[[nom_variable]],n.dens = n.dens,
                              nom=nom_variable,type_entree=type_entree,
                              dates_prises=dates_import,j_show = 19)
    l_AD[[nom_variable]]<-AN_GPD
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
  for(nom_variable in liste_noms){
    Donnes<-l_Orig_all[[nom_variable]]
    colnames(Donnes)<-c(1:37)
    rownames(Donnes)<-c(1:nrow(Donnes))
    Vecteur_DIM<-c(1:ncol(Donnes))
    # 1) Use GAM to analyse short-tail distribs ----------------------------------
    Donnes_pos<-Donnes[IND_select,]
    l_Orig[[nom_variable]]<-Donnes_pos
    Unif_EGPD<-matrix(NA,nrow = length(IND_select),
                      ncol=ncol(Donnes))
    print(nom_variable)
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
    # Other methods with transf and then heavy tail with Naveau code.

    # PARETO<-do.call(cbind.data.frame,Resultat_P$obs)
    # colnames(PARETO)<-c(1:ncol(PARETO))
    # CHOSEN_XI<-(1/2)
    # TRANSFO<-log(PARETO**CHOSEN_XI)
    # SIGMA<-1
    # INIT_GPD<-c(SIGMA,CHOSEN_XI)
    # for(t in c(1:length(Resultat_P$gamma))){
    #   VAR_t<-TRANSFO[,t]
    #   Model_EXTGP<-mev::fit.extgp(data =VAR_t,
    #                               model = 1,
    #                               method = "mle",
    #                               init =c(0.5,INIT_GPD),
    #                               plots = FALSE)
    #   Theta_extgp_t<-Model_EXTGP$fit$mle
    #   
    #   # Take estimation ---------------------------------------------------------
    #   UNIF_t<-sapply(X =VAR_t,
    #                  FUN=function(x){
    #                    mev::pextgp(q =x ,kappa =Theta_extgp_t[["kappa"]],
    #                                xi =Theta_extgp_t[["xi"]] ,
    #                                sigma = Theta_extgp_t[["sigma"]],
    #                                type = 1)})
    #   MATRIX_TRANSFO[,t]<-UNIF_t
    #   
    # }
    if(opt_Frech){
      FINAL_transfo<--1/log(Unif_EGPD)
    }else{
      FINAL_transfo<-1/(1-Unif_EGPD)
    }
    l_ALL[[nom_variable]]<-list("transf"=FINAL_transfo,
                                "params_transfo"=Vect_egpd,
                                "Fitting"=Vect_fitting,
                                "INIT"=Vect_Init)
    df[[nom_variable]]<-apply(X = FINAL_transfo,MARGIN = 1,
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
                             f_transf_inv){
  Array_simul<-array(data = NA,dim = c(length(l_variables),M,37))
  
  LISTE_obs_exts<-list()
  Liste_all<-result_transformation$resume
  Vect_l_function<-result_transformation$df
  
  # CF Kokozka ---------------------------------------------------------------
  L<-300
  vect_k<-c(20:L)
  l_prime<-apply(Vect_l_function,MARGIN = 1,
                 FUN = List_Params_RF[["function"]])
  Graphics_estimators_gamma(series = l_prime,
                            vect_k = vect_k,
                            Title_graphic = "Shape parameter of l_g",
                            NB_years = 37)
  
  # Choice of the individuals -----------------------------------------------
  Seuil_lprime<-quantile(l_prime,Q_thresh)
  
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
  Indices_exts<-which(l_prime>Seuil_lprime)
  #Rank transfo
  N<-nrow(Vect_l_function)
  Vect_l_orig<-Vect_l_function
  Vect_l_function<-apply(X = Vect_l_function,MARGIN = 2,FUN = function(x){
    Denom<-N+1-rank(x)
    return(N/Denom)
  })
  Excedents_lprime<-Vect_l_function[Indices_exts,]
  
  Excesses_maximum<-l_prime[Indices_exts]
  # Anderson-Darling test for the maximum -- Pareto ?  ----------------------
  print("p value AD test Pareto")
  print(goftest::ad.test(x =Excesses_maximum, 
                         null = extRemes::"pevd",
                         shape=1,scale=Seuil_lprime,
                         threshold=Seuil_lprime,
                         type="GP")$p.value)
  Theta_opt<-NA
  if(Params_risk_Function[["parametric"]]==TRUE){
    Theta_opt<-Estim_param_RF_homogeneous(Params_risk_Function= Params_risk_Function,
                                          Seuil_lprime= Seuil_lprime,
                                          Vect_l_function = Vect_l_function)
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
  Matrix_C<-function_Structure_Matrice(NB_dim = NbScores_Omega)
  Modele_coords<-VineCopula::RVineCopSelect(Unif_coord,
                                            Matrix = Matrix_C)
  Bool_filled<-FALSE
  N_sim<-M
  Coords_simul<-matrix(NA,nrow = M,ncol = NbScores_Omega)
  while(Bool_filled==FALSE){
    # Simulations of coordinates  ------------------------------------------------------------
    Simulations_coord<-VineCopula::RVineSim(N=N_sim,
                                            RVM = Modele_coords)
    colnames(Simulations_coord)<-1:ncol(Simulations_coord)
    Coords_ech_orig<-matrix(NA,nrow=N_sim,
                            ncol=NbScores_Omega)
    for(j in c(1:ncol(Simulations_coord))){
      sortie<-quantile(Scores[,j],Simulations_coord[,j])
      Coords_ech_orig[,j]<-as.numeric(sortie)
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
    if(Params_risk_Function[["parametric"]]==FALSE){
      Sim_l<-Simul_RF_max_NP(Params_risk_Function = Params_risk_Function,
                             M=N_sim,Seuil_lprime = Seuil_lprime,
                             Vect_l_function = Vect_l_function)
    }else{ 
      
      Sim_l<-Simulation_gParetoP(Params_risk_Function=Params_risk_Function,
                                 Threshold= Seuil_lprime,d=ncol(Vect_l_function),
                                 M=N_sim,theta_opt=Theta_opt)
      ### From Pareto to original.
      TAU<-0.95
      Sim_l<-sapply(X = c(1:ncol(Sim_l)),
                    FUN = function(x){
                      result<-rep(NA,nrow(Sim_l))
                      col_j<-Sim_l[,x]
                      data_j<-Vect_l_orig[,x]
                      Qj<-quantile(data_j,TAU)
                      ### MEV produces vectors with Frechet margins.
                      Unif_<-exp(-(col_j)^(-1))
                      Inds_exts<-which(Unif_>TAU)
                      Convert_todata<-Unif_
                      Sub_non_exts<-Unif_[-Inds_exts]
                      #return the corresponding data quantile for non extremes.
                      Convert_todata[-Inds_exts]<-sapply(X = Sub_non_exts,
                                      function(x){return(quantile(data_j,x))})
                      # And use the Pareto tails for the extreme ones.
                      V<-(Unif_-TAU)/(1-TAU)
                      Convert_todata[Inds_exts]<-Seuil_lprime*(1-V)^(-1)
                      return(Convert_todata)
                    })
      return(list("simuls"=Sim_l,"obs"=Vect_l_orig[Indices_exts,],
                  "th"=Seuil_lprime))
    }
    Sim_l<-as.data.frame(Sim_l)
    colnames(Sim_l)<-l_variables
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
                         apply(X = Z_varj,MARGIN = 1,FUN = calcul_norme_L2))
      indicatrice_pos<-which(apply(X=Z_varj,
                                   FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
      Unif<-exp(-Z_varj[indicatrice_pos,]^(-1))
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
      UNIF_k<-exp(-Z_vark^(-1))
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
      ### Unif-->Frechet to compare in Frechet scale of obs.
      LISTE_Frechet_SIMUL[[name_variable_for_conv]]<--1/log(Simul_whole_kunif)
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
      INDICES_ACP<-1:nrow(Simul_whole_kunif)
      ### Unif-->Frechet to compare in Frechet scale of obs.
      LISTE_Frechet_SIMUL[[name_variable_for_conv]]<--1/log(Simul_whole_kunif)
      NO_ACP<-lapply(INDICES_ACP,FUN = fnct_select_colonne,
                     df=Simul_whole_kunif)
      
      # Reconversion in the good scale--------------------------------------
      Theta_EXTGPD_k<-Liste_all[[name_variable_for_conv]]$params_transfo
      Variables_reconversion_ACP<-lapply(NO_ACP,function_reconversion_Unif_EXTGPD,
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
  return(list("Risk_MV"=l_prime[Indices_exts],
              "Param_found"=Theta_opt,"simul"=LISTE_simul,
              "obs_exts"=LISTE_obs_exts,"Indices_exts"=Indices_exts,
              "coords_simul"=Coords_simul, "coords_data"=Scores,
              "Frechet_normal"=list("obs"=LISTE_Frechet_OBS,
                                    "simul"=LISTE_Frechet_SIMUL),
              "EIGEN_functions"= EIGEN_functions,"family_copula"=Modele_coords,
              "threshold"=Seuil_lprime))
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

Inner_k_l<-function(l_kl,Mat){
  k<-l_kl[1]
  l<-l_kl[2]
  Vk<-Mat[,k]
  Vl<-Mat[,l]
  Prod_scal<-as.numeric(t(Vk)%*%Vk/length(Vk))
  return(Prod_scal)
}

Extreme_corr<-function(j,liste_MV_simul,l_name,L,d){
  Mat<-matrix(NA,nrow = L,
              ncol = d)
  for(l in c(1:d)){
    Name<-l_name[l]
    Df<-liste_MV_simul[[Name]]
   # Df<-sapply(Df,unlist)
    Mat[,l]<-unlist(Df[j,])
  }
  vect_pair<-list()
  i<-1
  for(l in c(1:d)){
    rest_<-c(1:d)
    one_<-rest_[l]
    rest_<-rest_[-l]
    for(l_tilde in c(1:length(rest_))){
      other<-rest_[l_tilde]
      pair<-c(one_,other)
      Reponse_bool<-sapply(vect_pair,function(x,ref){
        return(all(ref%in%x))
      },ref=pair)
      cond<-any(TRUE%in%Reponse_bool)
      if(cond==FALSE){
        vect_pair[[i]]<-pair
        i<-i+1
      }
    }
  }
  Num<-sapply(vect_pair,Inner_k_l,Mat=Mat)
  # Frac --------------------------------------------------------------------
  Denom<-max(apply(X = Mat,MARGIN = 2,FUN = calcul_norme_L2))^(2)
  return(Num/Denom)
}

Conv_scalar_product<-function(list_TS,Nb_inds_exts,l_name_variables,
                              g_function){
  L<-length(l_name_variables)
  L2_d1<-apply(list_TS[[l_name_variables[1]]],MARGIN = 1,
               FUN = calcul_norme_L2)
  M<-length(L2_d1)
  Base_l<-matrix(NA,nrow = M,
                 ncol=L)
  Base_l[,1]<-L2_d1
  for(j in c(2:L)){
    L2_dj<-apply(list_TS[[l_name_variables[j]]],MARGIN = 1,
                   FUN = calcul_norme_L2)
    Base_l[,j]<-L2_dj
  }
  Lg<-apply(X = Base_l,
               MARGIN = 1,
               FUN = g_function)
  M<-length(L2_dj)
  Qlevel<-1-(Nb_inds_exts/M)
  Inds_extremes<-which(Lg>quantile(Lg,Qlevel))
  list_extj<-list()
  list_conv_j<-list()
  for(j in c(1:L)){
    L2_extj<-Base_l[Inds_extremes,j]
    Ext_j<-list_TS[[l_name_variables[j]]][Inds_extremes,]
    #convergence of first moments.
    Angle_extj<-t(t(Ext_j)
                  %*%diag(L2_extj^(-1)))
    list_extj[[l_name_variables[j]]]<-Ext_j
    list_conv_j[[l_name_variables[j]]]<-Min_function_conv(
      Shape_d =Angle_extj )
  }
  #convergence of scalar products.
  Sigmaxy<-sapply(c(1:length(Inds_extremes)),
               liste_MV_simul = list_extj,
               Extreme_corr,
               l_name =l_name_variables,
               L = ncol(Ext_j),
               d = length(l_name_variables))
  return(list("corr"=mean(Sigmaxy),
              "sd_corr"=sd(Sigmaxy),
              "first_moments"=list_conv_j))
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

