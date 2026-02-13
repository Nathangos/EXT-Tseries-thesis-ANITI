Apply_VAR_per_t<-function(t,Epsi_sim,X_delta,Indexes_delta,ARRAY_VAR){
  MOD_VAR_t<-ARRAY_VAR[t,,]
  nb_vars<-ncol(MOD_VAR_t)
  MAT_beta<-MOD_VAR_t[c(1:nb_vars),c(1:nb_vars)]
  X_delta_t<-X_delta[Indexes_delta,t,]
  AX<-X_delta_t%*%(MAT_beta)
  const_index<-nb_vars+1
  Coeff_tend_<-MOD_VAR_t[const_index,]
  Const_<-MOD_VAR_t[const_index,]
  Rep_CONST<-matrix(rep(Const_,nrow(AX)),
                    nrow =nrow(AX),ncol(AX),
                    byrow = TRUE) 
  Coeff_tend_<-MOD_VAR_t[nrow(MOD_VAR_t),]
  TEND_<-Indexes_delta*Coeff_tend_
  Determ_part<-AX+TEND_+Rep_CONST
  Random_part<-Epsi_sim[,t, ]
  Predictions<-Determ_part+Random_part
  return(Predictions)
}

# Knn with choice of the index --------------------------------------------

KNN_distancef<-function(x,X_data,f_distance,K){
  DL2<-apply(X=X_data,FUN=function(x_ref,x){
    # Relative Squared error to the simulated since the scale are not the same
    return(sum((x-x_ref/x_ref)**2))},x_ref=x,
    MARGIN = 1)
  Indexes_K_closest<-order(DL2)[1:K]
  return(sample(x = Indexes_K_closest,size = 1,
                prob = rep(1/K,K)))
}

Generator_TSERIES_MV<-function(coeurs,type_donnees,list_variable,
                               Name_riskF,link_import_residuals,
                               link_VAR_model,
                               K_chosen,f_distance,link_export_generations,
                               cols_gg){
  link_g_chosen<-paste0("residuals_MV/",Name_riskF,"/")
  # Import X(M-1) -----------------------------------------------------------
  ######
  repertory<-"../ss_tend/"
  if(type_donnees=="HIVER"){
    repertory<-paste0(repertory,"HIVER/")
  }
  l_Orig<-list()
  l_Epsi_data<-list()
  l_Transf<-list()
  l_RiskF<-list()
  l_Epsi_Sim<-list()
  l_RiskF_Epsi_Data<-list()
  l_RiskF_Epsi_Sim<-list()
  TS_MV_Predictions<-list()
  d<-length(list_variable)
  l_Epsi_DOrig<-list()
  for(name_variable in list_variable){
    lien_donnees<-paste0(repertory,name_variable,"_ss_tend.csv")
    Donnes<-read.csv(file=lien_donnees)[,2:38]
    colnames(Donnes)<-c(1:37)
    rownames(Donnes)<-c(1:nrow(Donnes))
    l_Orig[[name_variable]]<-Donnes
    N_data<-nrow(Donnes)
    C_data<-ncol(Donnes)
    l_RiskF[[name_variable]]<-apply(X = Donnes,MARGIN = 1,
                                    FUN = calcul_norm_L2)
    # 1) Import Resid data -------------------------------------------------------
    Epsi_data<-read.csv(file=paste0(link_import_residuals,
                                    name_variable,"_residuals.csv"))[,c(2:38)]
    colnames(Epsi_data)<-c(1:ncol(Epsi_data))
    l_Epsi_DOrig[[name_variable]]<-Epsi_data
    # Drop the first resid
    Epsi_data<-Epsi_data[c(2:nrow(Epsi_data)),]
    l_Epsi_data[[name_variable]]<-Epsi_data
    l_RiskF_Epsi_Data[[name_variable]]<-apply(X = Epsi_data,
                                              MARGIN = 1,
                                              FUN=calcul_norm_L2)
    
    # 2)Import Resid simulated ------------------------------------------------------------
    Epsi_Sim<-read.csv(file = paste0(link_g_chosen,
                                     name_variable,"_simul.csv"))[,c(2:38)]
    colnames(Epsi_Sim)<-colnames(Epsi_data)
    N_sim<-nrow(Epsi_Sim)
    l_Epsi_Sim[[name_variable]]<-Epsi_Sim
    l_RiskF_Epsi_Sim[[name_variable]]<-apply(X = Epsi_Sim,
                                             MARGIN = 1,
                                             FUN = calcul_norm_L2)
  }
  # 2) Import indexes of exceedances -------------------------------------------
  Indexes_exceed<-read.csv(file=paste0(link_g_chosen,
                                       "ext_series_chosen_risks.csv"))
  Precise_index<-Indexes_exceed[,3]
  Matrix_Epsi_sim<-array(NA,dim = c(N_sim,C_data,d))
  Matrix_X_delta<-array(NA,dim = c(N_data,C_data,d))
  Matrix_EXTREMES_events<-array(NA,dim = c(nrow(Indexes_exceed),
                                           C_data,d))
  for(Z in c(1:d)){
    name_variable<-list_variable[Z]
    ELT<-as.matrix(l_Orig[[name_variable]])
    ELT_Epsi_sim<-as.matrix(l_Epsi_Sim[[name_variable]])
    Matrix_Epsi_sim[,,Z]<-ELT_Epsi_sim
    Matrix_X_delta[,,Z]<-ELT
    Matrix_EXTREMES_events[,,Z]<-ELT[Precise_index,]
  }
  RiskF<-as.data.frame(l_RiskF)
  RiskF_Epsi_Sim<-as.data.frame(l_RiskF_Epsi_Sim)
  RiskF_Epsi_Data<-as.data.frame(l_RiskF_Epsi_Data)
 
  # 3) Import VAR model
  END_column<-length(list_variable)+1
  VAR_model<-read.csv(file = link_VAR_model)[,c(2:END_column)]
  VAR_model<-as.data.frame(VAR_model)
  ### Prevent order problems.
  VAR_model_adjusted<-list()
  for(name_var in lNAME){
    VAR_model_adjusted[[name_var]]<-VAR_model[[name_var]]
  }
  VAR_model_adjusted<-as.data.frame(VAR_model_adjusted)
  # 4) Choose states
  Indexes_knn_sim_vs_data<-apply(X = RiskF_Epsi_Sim,MARGIN = 1,
                                 FUN = KNN_distancef,f_distance=f_distance,
                                 K=K_chosen,X_data=RiskF_Epsi_Data)
  Index_delta<-Indexes_knn_sim_vs_data-1
  L<-nrow(RiskF_Epsi_Data)
  END<-L-1
  Z_gg<-1
  list_gg<-list()
  list_gg1<-list()
  for(name_v in list_variable){
    Data_Xo_v<-l_RiskF[[name_v]][1:END]
    Epsi_data_v<-l_RiskF_Epsi_Data[[name_v]]
    Epsi_Sim_v<-l_RiskF_Epsi_Sim[[name_v]]
    Data_Xo_chosen<-Data_Xo_v[Index_delta]
    Epsi_data_v<-Epsi_data_v[2:L]
    MATRIX_rel_epsi_previous<-cbind.data.frame(Data_Xo_v,Epsi_data_v)
    colnames(MATRIX_rel_epsi_previous)<-c("X_prevnorm","Epsi_norm")
    MATRIX_rel_epsi_previous_sim<-cbind.data.frame(Data_Xo_chosen,Epsi_Sim_v)
    colnames(MATRIX_rel_epsi_previous_sim)<-c("X_prevnorm","Epsi_norm")
    m<-min(apply(X = MATRIX_rel_epsi_previous_sim,MARGIN = 2,FUN = min),
           apply(X = MATRIX_rel_epsi_previous,MARGIN = 2,FUN = min))
    M<-max(apply(X = MATRIX_rel_epsi_previous_sim,MARGIN = 2,FUN = max),
           apply(X = MATRIX_rel_epsi_previous,MARGIN = 2,FUN = max))
    XLAB<-expression("L2 norm of "~epsilon[M~","~Z_gg])
    XLAB<-as.expression(do.call('substitute', list( XLAB[[1]], list(Z_gg=Z_gg))))
    YLAB<-expression("L2 norm of "~tilde(X)[M-Delta~","~Z_gg])
    YLAB<-as.expression(do.call('substitute', list( YLAB[[1]], list(Z_gg=Z_gg))))
    GG_relation<-ggplot(MATRIX_rel_epsi_previous,aes(x=Epsi_norm,y=X_prevnorm,
                                                     col="data"))+
      geom_point()+
      geom_point(data=MATRIX_rel_epsi_previous_sim,aes(y=X_prevnorm,
                 x=Epsi_norm,col="simulations",shape="simulations"),size=1
                 ,pch=17)+
      xlab(XLAB)+
      ylab(YLAB)+
      scale_color_manual(values=cols_gg)+
      ylim(c(m,M))+
      xlim(c(m,M))+
      labs(col="Legend")
    list_gg[[name_v]]<-GG_relation
    Z_gg<-Z_gg+1

  }
  rel_MV_epsi_previous<-do.call("grid.arrange", c(list_gg,
                                                  ncol=length(list_variable)))
  ggsave(plot= rel_MV_epsi_previous,width=8,height=6,
         filename=paste0("graphiques_MV/data/",Name_riskF,"/rel_Epsi_previous_observation.png"))
  
  # 5) Apply VAR model
  # Matrix to Array
  p<-nrow(VAR_model)/C_data
  ARRAY_VARt<-array(NA,dim=c(C_data,p,d))
  for(T in c(1:C_data)){
    BEG<-(T-1)*p+1
    END<-T*p
    ELT<-VAR_model_adjusted[c(BEG:END),]
    ARRAY_VARt[T,,]<-as.matrix(ELT)
  }
  # Function by time
  Vect_times<-c(1:C_data)
  X_pred_per_time<-lapply(X = Vect_times,FUN=Apply_VAR_per_t,Epsi_sim=Matrix_Epsi_sim,
                          X_delta=Matrix_X_delta,Indexes_delta=Index_delta,
                          ARRAY_VAR=ARRAY_VARt)
  
  # 6) Reorganise according to variable for future export
  Matrix_X_end<-array(NA,dim = c(d,N_sim,C_data))
  Matrix_Previous<-array(NA,dim = c(d,N_sim,C_data))
  for(T in c(1:C_data)){
    for(k in c(1:d)){
      Matrix_X_end[k,,T]<-X_pred_per_time[[T]][,k]
    }
    
  }
  # previous obs
  for(Zexport in c(1:d)){
    name_variable<-list_variable[Zexport]
    ELT<-as.matrix(l_Orig[[name_variable]])
    Matrix_Previous[Zexport,,]<-ELT[Index_delta,]
  }
  # Export per variable
  LIST_DELTA<-list()
  LIST_TARGET<-list()
  LIST_Result<-list()
  for(Zexport in c(1:d)){
    name_forcing<-lNAME[Zexport]
    SUB<-Matrix_X_end[Zexport,,]
    LIST_Result[[name_forcing]]<-SUB
    write.csv(x=SUB,
              file = paste0(link_export_generations,
                            name_forcing,"_fsim.csv"))
    ## Export the chosen time series
    ELT_previous_k<-Matrix_Previous[Zexport,,]
    LIST_DELTA[[name_forcing]]<-ELT_previous_k
    write.csv(x=ELT_previous_k,
              file = paste0(link_export_generations,
                            name_forcing,"_f_previous.csv"))
    
    ## Export extreme events
    SUB_data<-l_Epsi_DOrig[[name_forcing]]
    SUB_data<-SUB_data[Precise_index,]
    LIST_TARGET[[name_forcing]]<-SUB_data
    write.csv(x= SUB_data,
              file = paste0(link_export_generations,
                            name_forcing,"_fdata.csv"))
    #Matrix_EXTREMES_events[,,Z]<-ELT[Precise_index,]
  }
  return(list("X0"=LIST_DELTA,
              "Epsi_sim"=l_Epsi_Sim,
              "Target"=LIST_TARGET,
              "Result"=LIST_Result,
              "Epsi_Data"=l_Epsi_data))
}
