# Apply for each time the corresponding VAR model--------------------------------------------
#####################

#' Apply_VAR_per_t
#'
#' @param t int. Time step of the multivariate time series
#' @param list_sim list[str: vector]. For each time 
#' step, the simulated residual vector of variable values.
#' @param listX_delta  list[str: vector]. For each time 
#' step, the vector of previous multivariate time series X(M-Delta)
#' @param list_VAR list[str]. Each value contains VAR parameters.
#' @param Indexes_delta vector[int]. Indexes of chosen previous
#' multivariate time series
#'
#' @return
#' @export
#'
#' @examples
Apply_VAR_per_t<-function(t,list_sim,listX_delta,
                          list_VAR,Indexes_delta){
  
  ### Pass by index or value
  ### index --> choice previous returning index
  ### value --> choice previous returning a fictive time series (KNN_2)
  t_str<-as.character(t)
  MOD_VAR_t<-as.data.frame(list_VAR[[t]])
  Var_dimension<-ncol(MOD_VAR_t)
  
  END_reg<-Var_dimension-1
  Const_<-as.numeric(MOD_VAR_t[Var_dimension,c(2:Var_dimension)])
  Beta_params<-MOD_VAR_t[c(1:END_reg),c(2:Var_dimension)]
  MAT_beta<-as.matrix(Beta_params)
  Order_variables<-colnames(Beta_params)
  
  X_delta_t_whole<-as.data.frame(listX_delta[[t_str]])
  if(is.null(Indexes_delta)){
    X_delta_t<-as.matrix(X_delta_t_whole[,Order_variables])
  }else{
    X_delta_t<-as.matrix(X_delta_t_whole[Indexes_delta,
                                         Order_variables])
  }
  ### Use the same ordering for the previous observation
  ### here --> if we use indexes delta, from =l_X_data
  ### otherwise--> use directly another argument--> NEW_list[[t_str]]
  
  ###Regression
  AX<-X_delta_t%*%MAT_beta
  ### Constant part
  Rep_CONST<-matrix(rep(Const_,nrow(AX)),
                    nrow =nrow(AX),ncol(AX),
                    byrow = TRUE) 
  Determ_part<-AX+Rep_CONST

  ## Epsilon simulated to add
  Random_part<-as.data.frame(list_sim[[t_str]])
  Random_part<-as.matrix(Random_part[,Order_variables])
  Predictions<-Determ_part+Random_part
  rownames(X_delta_t)<-c(1:nrow(X_delta_t))
  rownames(Predictions)<-c(1:nrow(Predictions))
  return(list("simulated"=Predictions,
              "previous"=X_delta_t))
}

distance_sq<-function(h){
  return(sum((h)^2)**(1/2))
}
# Knn with choice of the index --------------------------------------------
#####################
KNN_distancef<-function(i,X_sim,X_data,f_distance,K,hearts){
  
  ### Compute distance for each forcing condition
  if(is.null(f_distance)){
  
    ### freely choosing the previous observation
    Nber_candidates<-nrow(X_data[[1]])
    Sampled_M<-sample(Nber_candidates,
                      size = 1)
  }else{
    ### using a criterion to choose it.
    if(!is.null(names(X_data))){
      Variables_seen<-names(X_data)
      ### if we provide a matrix
      Varwise_dces<-sapply(Variables_seen,function(x){
        Varj_data<-as.matrix(X_data[[x]])
        Varj_sim<-X_sim[[x]]
        Varij_sim<-unlist(Varj_sim[i,])
        Gap_mat<-sweep(Varj_data, 2, Varij_sim, "-")
        Dces_found<-apply(Gap_mat,MARGIN = 1,
                          FUN = f_distance)
        return(Dces_found)
      })
      Summary_varwise_dces<-rowMeans(Varwise_dces)
    }else{
      ### if we provide a vector of gl
      D<-X_data-X_sim[i]
      Summary_varwise_dces<-sapply(D,FUN = f_distance)
    }
    
    ### Find closest time series
    sorted_dce<-sort(Summary_varwise_dces)[1:K]
    Indexes_K_closest<-order(Summary_varwise_dces)[1:K]
    # weight_per_k<-1/sorted_dce
    # W<-weight_per_k/sum(weight_per_k)
    weight_per_k<-1/sorted_dce
    W<-weight_per_k/sum(weight_per_k)
    Sampled_M<-sample(x = Indexes_K_closest,size = 1,
                      prob = W)
  }
  
  return(Sampled_M)
}
KNN_distancef_returnMean<-function(i,X_sim,X_data,f_distance,K,hearts,
                                   Y_data){
  Variables_seen<-names(X_sim)
  ### Compute distance for each forcing condition
  if(is.null(f_distance)){
    ### freely choosing the previous observation
    Nber_candidates<-nrow(X_data[[1]])
    throw("Invalid argument")
  }else{
    ### using a criterion to choose it.
    Varwise_dces<-sapply(Variables_seen,function(x){
      Varj_data<-as.matrix(X_data[[x]])
      ### Take into account the scale of the forcing condition
      Weight_summarised<-sum(Varj_data**2)**(-1/2)
      Varj_sim<-X_sim[[x]]
      if(is.vector(Varj_sim)){
        Varij_sim<-Varj_sim[i]
      }else{
        Varij_sim<-Varj_sim[i,]
      }
      Varij_sim<-unlist(Varij_sim)
      Gap_mat<-sweep(Varj_data, 2, Varij_sim, "-")
      Dces_found<-apply(Gap_mat,MARGIN = 1,
                        FUN = f_distance)
      #*Weight_summarised
      return(Dces_found)
    })
    Summary_varwise_dces<-rowMeans(Varwise_dces)
    ### Find closest time series
    sorted_dce<-sort(Summary_varwise_dces)[1:K]
    Indexes_K_closest<-order(Summary_varwise_dces)[1:K]
    # weight_per_k<-1/sorted_dce
    # W<-weight_per_k/sum(weight_per_k)
    # considering the distance
    # if(sum(sorted_dce)!=0){
    #   weight_per_k<-(sorted_dce)^(-1)
    #   W<-weight_per_k/sum(weight_per_k)
    # }else{
    #   W<-rep(1/length(Indexes_K_closest),length(Indexes_K_closest))
    # }
    # without considering the distance
    W<-rep(1/length(Indexes_K_closest),length(Indexes_K_closest))
    print(W)
    Output<-lapply(Variables_seen,function(x){
      ### take the neighbors' output
      Varj_Ydata<-as.matrix(Y_data[[x]])
      Selected_ones<-Varj_Ydata[Indexes_K_closest,]
      SUM_W<-t(Selected_ones)%*%W
      return(t(SUM_W)[1,])
    })
    names(Output)<-Variables_seen
    return(Output)
  }

}
Generator_TSERIES_MV<-function(hearts,list_variable,
                               Name_riskF,link_import_residuals,
                               link_rootVAR_model, K_chosen,
                               f_distance,link_export_generations,
                               cols_gg,prefix_link_resid,
                               subfix_link_Sim,Dates_Johanna,
                               subfix_link_Data,Vect_VARtimes,
                               Risk_Function,opt_used,
                               repertory_data){
  
  link_export_SIM<-paste0(link_export_generations,opt_used,"_")
  link_g_chosen<-paste0("residuals_MV/",Name_riskF,"/")

  # (1) Import X(M-1) -----------------------------------------------------------
  ######
  l_RiskF<-list()
  l_Orig<-lapply(Vect_VARtimes,
                    function(x){
                      return(list())
                    })
  l_Transf<-list()
  l_X_data<-list()
  l_Epsi_Sim<-l_X_data
  l_Epsi_Sim_for_export<-list()
  l_Epsi_Prev_Resid_for_export<-list()
  
  l_Epsi_DOrig<-list()
  l_RiskF_Epsi_Data<-list()
  l_RiskF_Epsi_Sim<-list()
  TS_MV_Predictions<-list()
  d<-length(list_variable)
  l_Epsi_DOrig<-list()
  
  l_Epsi_Sim_varwise<-list()
  l_Epsi_RSim_varwise<-list()
  
  l_Epsi_Data_varwise<-list()
  l_Epsi_RData_varwise<-list()
  
  # (2) Import indexes of exceedances -------------------------------------------
  #########
  Indexes_exceed<-read.csv(file=paste0(link_g_chosen,
                                       "ext_series_chosen_risks_",
                                       d,".csv"))
  Precise_index<-Indexes_exceed[,3]+1
  l_Epsi<-list()
  l_Target_Joh<-list()
  Params_scale<-list()
  Params_mu<-list()
  for(name_variable in list_variable){
    if(name_variable=="U"){
      link_transit<-paste0(repertory_data,name_variable)
      link_data<-paste0(link_transit,"_trunc_detrend.csv")
      Data<-read.csv(file=link_data)
      
    }else{
      link_transit<-paste0(repertory_data,name_variable)
      link_data<-paste0(link_transit,"_detrend.csv")
      Data<-read.csv(file=link_data)
    }
    D<-ncol(Data)
    Data<-Data[,c(2:D)]
    colnames(Data)<-substr(x = colnames(Data),
                             start = 2,stop = 3)
    rownames(Data)<-c(1:nrow(Data))
    l_Orig[[name_variable]]<-Data
    ind_target_joh<-Dates_Johanna+1
    l_Target_Joh[[name_variable]]<-Data[ind_target_joh,]
    N_data<-nrow(Data)
    C_data<-ncol(Data)
    Max_X_raw<-apply(X = Data,MARGIN = 1,
                     FUN=max)
    l_RiskF[[name_variable]]<-Max_X_raw
    
    # Import Resid data -------------------------------------------------------
    Epsi_data<-read.csv(file=paste0(link_import_residuals,
                                    name_variable,"_residuals_",
                                d,".csv"))
    colnames(Epsi_data)<-substr(x = colnames(Epsi_data),
                                start = 2,stop = 3)
    Mu_epsi_data<-colMeans(Epsi_data)
    Sd_epsi_data<-apply(Epsi_data,MARGIN = 2,FUN = sd)
    ### Used for the window
    L2_epsi_data<-apply(Epsi_data,MARGIN = 1,
                      FUN =calcul_norm_L2)
    Mu_epsi_Rdata<-mean(L2_epsi_data)
    Sd_epsi_Rdata<-sd(L2_epsi_data)
    l_Epsi_Data_varwise[[name_variable]]<-scale(Epsi_data)
    
    ### For window
    l_Epsi_RData_varwise[[name_variable]]<-scale(L2_epsi_data)
    
    Epsi_EXT<-Epsi_data[Precise_index,]
    l_Epsi_Prev_Resid_for_export[[name_variable]]<-Epsi_EXT
    
    # Drop the first resid if it was done before...
    ####
    if(nrow(Epsi_data)==nrow(Data)){
      Epsi_data<-Epsi_data[c(2:nrow(Epsi_data)),]
    }
    Max_epsiraw<-apply(X = Epsi_data,
                   MARGIN = 1,
                   FUN=max)
    l_RiskF_Epsi_Data[[name_variable]]<-Max_epsiraw
    # Import Resid simulated ------------------------------------------------------------
    Epsi_Sim<-read.csv(file = paste0(link_g_chosen,prefix_link_resid,
                                     name_variable,"_simul_",
                                     subfix_link_Sim,".csv"))
    colnames(Epsi_Sim)<-colnames(Epsi_data)
    Max_EpsiSim<-apply(X = Epsi_Sim,
                   MARGIN = 1,
                   FUN = max)
    l_Epsi_Sim_for_export[[name_variable]]<-Epsi_Sim
    l_RiskF_Epsi_Sim[[name_variable]]<-Max_EpsiSim
    l_Epsi_Sim_varwise[[name_variable]]<-scale(Epsi_Sim,
                               center =  Mu_epsi_data,
                               scale = Sd_epsi_data)
    
    ### For window
    L2_sim<-apply(X = Epsi_Sim,MARGIN = 1,FUN = calcul_norm_L2)
    l_Epsi_RSim_varwise[[name_variable]]<-scale(L2_sim,
                               center = Mu_epsi_Rdata,
                               scale = Sd_epsi_Rdata)
    
    Params_scale[[name_variable]]<-Sd_epsi_data
    Params_mu[[name_variable]]<-Mu_epsi_data
    for(Time in colnames(Epsi_Sim)){
      Add_data<-Data[,Time]
      Add_sim<-Epsi_Sim[,Time]
      l_Epsi_Sim[[Time]][[name_variable]]<-Add_sim
      l_X_data[[Time]][[name_variable]]<-Add_data
    }
  }
  
  # (3) Import the VAR models
  list_var_model<-list()
  for(elt in Vect_VARtimes){
    list_var_model[[elt]]<-read.csv(
      file = paste0(link_rootVAR_model,"VAR_model_",
                    elt,".csv"))
  }
  
  # (4) Choose states
  #l_RiskF_Epsi_Sim
  #l_RiskF_Epsi_Data
  Nsim<-nrow(l_Epsi_Sim_varwise[[1]])
  ### Input of KNN
  Df_RData<-do.call("cbind.data.frame",
                   l_Epsi_RData_varwise)
  Df_RSim<-do.call("cbind.data.frame",
                  l_Epsi_RSim_varwise)
  ### Compute the value of the gl in 
  ### simulations and in observations
  RF_data<-apply(X = Df_RData,MARGIN = 1,FUN = Risk_Function)
  RF_sim<-apply(X = Df_RSim,MARGIN = 1,FUN = Risk_Function)
  
  #Johanna
  Index_Joh_rep<-rep(Dates_Johanna,Nsim)
  if(opt_used!="KNN_2"){
    l_X_PREV<-l_X_data
    if(opt_used=="window"){
      print("window approach")
      Index_delta<-sapply(X =RF_sim,Sample_window,
                          y=NULL,x_window=RF_data,
                          size_window=20,use_weight=TRUE)
    }
    if(opt_used=="KNN"){
      print("KNN approach")
      #l_Epsi_RSim_varwise,l_Epsi_RData_varwise,l_Epsi_Sim_varwise,l_Epsi_Data_varwise
      Indexes_knn_sim_vs_data<-parallel::parSapply(cl = hearts,X =c(1:Nsim),
                           FUN = KNN_distancef,f_distance=f_distance,
                           K=K_chosen,X_data=l_Epsi_RData_varwise,
                           X_sim=l_Epsi_RSim_varwise,
                           hearts=hearts)
      Index_delta<-Indexes_knn_sim_vs_data
    }
    if(opt_used=="free"){
      print("Uncond approach")
      Index_delta<-sample(c(1:length(RF_data)),size = Nsim,
                          replace = TRUE)
    }
    #Johanna
    Index_Joh_rep<-rep(Dates_Johanna,length(Index_delta))
    #Index_delta<-Precise_index
    L<-nrow(l_Epsi_Data_varwise[[1]])
    END<-L-1
    Z_gg<-1
    list_gg1<-list()
    FSIZE<-15
    Corrected_list_variable<-sapply(list_variable,
                                    FUN = Fct_correct_name,
                                    target="Surcote",replacement="Surge")
    for(name_v in list_variable){
      Data_Xo_v<-l_RiskF[[name_v]][1:END]
      Epsi_data_v<-l_RiskF_Epsi_Data[[name_v]]
      Epsi_Sim_v<-l_RiskF_Epsi_Sim[[name_v]]
      ### if you don't use a real obs--> use directly output of KNN
      Data_Xo_chosen<-Data_Xo_v[Index_delta]
      Epsi_data_v<-Epsi_data_v[2:L]
      MATRIX_rel_epsi_previous<-cbind.data.frame(Data_Xo_v,
                                                 Epsi_data_v)
      colnames(MATRIX_rel_epsi_previous)<-c("X_prevnorm","Epsi_norm")
      MATRIX_rel_epsi_previous_sim<-cbind.data.frame(Data_Xo_chosen,
                                                     Epsi_Sim_v)
      colnames(MATRIX_rel_epsi_previous_sim)<-c("X_prevnorm","Epsi_norm")
      m<-min(apply(X = MATRIX_rel_epsi_previous_sim,MARGIN = 2,FUN = min),
             apply(X = MATRIX_rel_epsi_previous,MARGIN = 2,FUN = min))
      M<-max(apply(X = MATRIX_rel_epsi_previous_sim,MARGIN = 2,FUN = max),
             apply(X = MATRIX_rel_epsi_previous,MARGIN = 2,FUN = max))
      XLAB<-expression("Max of "~epsilon[M~","~Z_gg])
      XLAB<-as.expression(do.call('substitute', list( XLAB[[1]], 
                                                      list(Z_gg=Corrected_list_variable[Z_gg]))))
      YLAB<-expression("Max of "~tilde(X)[M-Delta~","~Z_gg])
      YLAB<-as.expression(do.call('substitute', list( YLAB[[1]], 
                                                      list(Z_gg=Corrected_list_variable[Z_gg]))))
      GG_relation<-ggplot2::ggplot(MATRIX_rel_epsi_previous,aes(x=Epsi_norm,y=X_prevnorm,
                                                       col="data"))+
        geom_point()+
        geom_point(data=MATRIX_rel_epsi_previous_sim,aes(y=X_prevnorm,
                                                         x=Epsi_norm,col="simulations",shape="simulations"),size=0.75
                   ,pch=17)+
        xlab(XLAB)+
        ylab(YLAB)+
        scale_color_manual(values=cols_gg)+
        # ylim(c(m,M))+
        # xlim(c(m,M))+
        labs(col="Legend")+
        theme(legend.direction = "horizontal",
              axis.title=element_text(size=FSIZE),
              strip.text=element_text(size=13),
              legend.title = element_text(size=14),
              legend.text=element_text(size=13))
      if(Z_gg==1){
        LGD_relation<-get_legend(GG_relation)
        Summary_relation<-GG_relation+
          guides(col="none")
      }else{
        Summary_relation<-Summary_relation+
          GG_relation+guides(col="none")
      }
      #list_gg[[name_v]]<-GG_relation
      Z_gg<-Z_gg+1
      
    }
    rel_MV_epsi_previous<-plot_grid(Summary_relation,LGD_relation,
                                    ncol=1,rel_heights =c(10,1))
    f_previous<-ggdraw()+
      draw_plot(rel_MV_epsi_previous)+
      theme(plot.background = element_rect(fill = "white",
                                           color = NA))
    ggsave(filename=paste0("Graphics_MV/data/",Name_riskF,
                           "/rel_Epsi_previous_observation.png"),
           plot = f_previous,width=10,height = 6)

  }else{
    ### Use a fictive time series (average of neighbors)--> do not use
    ### the index in the dataframe
    print("KNN2 approach")
    Index_delta<-NULL
    List_Delta_Avg<-parallel::parLapply(cl = hearts,X =c(1:Nsim),
                        fun = KNN_distancef_returnMean,f_distance=f_distance,
                        K=K_chosen,X_data=l_Epsi_Data_varwise,
                        X_sim=l_Epsi_Sim_varwise,hearts=hearts,
                        Y_data=l_Orig)
    Combined<-do.call("rbind",List_Delta_Avg)
    l_X_PREV<-list()
    for(number_col in c(1:ncol(Combined))){
      nameV_Fill<-colnames(Combined)[number_col]
      Df_v<-do.call("rbind",
                    Combined[,number_col])
      CNAMES<-colnames(Df_v)
      for(elt in CNAMES){
        Content_time_nv<-Df_v[,elt]
        l_X_PREV[[elt]][[nameV_Fill]]<-Content_time_nv
      }
    }
  }
  # (5) Apply VAR model
  #####################
  Vect_times<-c(1:C_data)
  X_simul_previous_per_time<-lapply(X = Vect_times,FUN=Apply_VAR_per_t,
                                    list_sim=l_Epsi_Sim,listX_delta=l_X_PREV,
                                    list_VAR=list_var_model,
                                    Indexes_delta=Index_delta)
 
  X_prev_Joh_per_time<-lapply(X = Vect_times,FUN=Apply_VAR_per_t,
                                    list_sim=l_Epsi_Sim,listX_delta=l_X_data,
                                    list_VAR=list_var_model,
                              Indexes_delta=Index_Joh_rep)
  names(X_simul_previous_per_time)<-as.character(Vect_times)
  names(X_prev_Joh_per_time)<-names(X_simul_previous_per_time)
  # (6) Reorganise according to variable for future export
  #################
  list_Xend<-list()
  list_XendJoh<-list()
  list_Xprevious<-list()
  for(T_seen in names(X_simul_previous_per_time)){
    Simulated_t<-X_simul_previous_per_time[[T_seen]][["simulated"]]
    Simulated_tJ<-X_prev_Joh_per_time[[T_seen]][["simulated"]]
    Previous_t<-X_simul_previous_per_time[[T_seen]][["previous"]]
    Var_used<-colnames(Simulated_t)
    for(Name_variable in Var_used){
      if(!Name_variable%in%names(list_Xend)){
        L<-list()
        L[[T_seen]]<-Simulated_t[,Name_variable]
        list_Xend[[Name_variable]]<-L
        
        ###Johanna
        L_Joh<-list()
        L_Joh[[T_seen]]<-Simulated_tJ[,Name_variable]
        list_XendJoh[[Name_variable]]<-L_Joh
        
        ###previous
        L_prev<-list()
        L_prev[[T_seen]]<-Previous_t[,Name_variable]
        list_Xprevious[[Name_variable]]<-L_prev
        
      }else{
        list_Xend[[Name_variable]][[T_seen]]<-Simulated_t[,Name_variable]
        list_Xprevious[[Name_variable]][[T_seen]]<-Previous_t[,Name_variable]
        
        ###Johanna
        list_XendJoh[[Name_variable]][[T_seen]]<-Simulated_tJ[,Name_variable]
      }
    }
    
    
  }
  list_Xend<-lapply(list_Xend,FUN = function(x){
    return(do.call(cbind.data.frame,x))
  })
  list_Xprevious<-lapply(list_Xprevious,FUN = function(x){
    return(do.call(cbind.data.frame,x))
  })
  list_XendJoh<-lapply(list_XendJoh,FUN = function(x){
    return(do.call(cbind.data.frame,x))
  })
  # Export per variable
  ######
  LIST_TARGET<-list()
  list_previous_Target<-list()
  for(Zexport in c(1:d)){
    ### Export for DoE/metamodeling steps the generated extremes
    ######
    name_forcing<-lNAME[Zexport]
    SUB<-list_Xend[[name_forcing]]
    write.csv(x=SUB,
              file = paste0(link_export_SIM,
                            prefix_link_resid,
                            name_forcing,"_fsim_",
                            subfix_link_Sim,".csv"),
              row.names = FALSE)
    ## Export the chosen time series
    Previous<-list_Xprevious[[name_forcing]]
    write.csv(x=Previous,
              file = paste0(link_export_generations,
                            name_forcing,"_f_previous_d=",
                            d,".csv"),row.names = FALSE)

    ## Export extreme events
    SUB_data<-l_Orig[[name_forcing]]
    SUB_data_ext<-SUB_data[Precise_index,]
    Precise_prev_index<-Precise_index-1
    SUB_data_prev_ext<-SUB_data[Precise_prev_index,]
    
    list_previous_Target[[name_forcing]]<-SUB_data_prev_ext
    LIST_TARGET[[name_forcing]]<-SUB_data_ext
    write.csv(x= SUB_data_ext,
              file = paste0(link_export_generations,
                            name_forcing,"_fdata_",
                            subfix_link_Data,".csv"),
              row.names = FALSE)
  }
  return(list("X0"=list_Xprevious,"Epsi_sim"=l_Epsi_Sim_for_export,
              "Epsi_data_ext"=l_Epsi_Prev_Resid_for_export,
              "Target"=LIST_TARGET,"previous_target"=
                list_previous_Target,"Result"=list_Xend,
              "index_used"=Index_delta,"XwholeObs"=l_X_data,
              "predictions_Johanna"=list_XendJoh,"target_Johanna"=l_Target_Joh))
  
}
K_fold_k_param_j_test<-function(hearts,f_distance,l_Epsi_Data_varwise,
                              l_output,indexes_Kfold,j_test,K_chosen){
  Indexes_train<-which(indexes_Kfold!=j_test)

  l_Epsi_Data_varwise_train<-list()
  l_Epsi_Data_varwise_test<-list()
  
  l_output_train<-list()
  l_output_test<-list()

  for(nameV in names(l_Epsi_Data_varwise)){
    ### Input of Knn
    ######
    Df_raw<-l_Epsi_Data_varwise[[nameV]]
    if(is.vector(Df_raw)){
      Df_train<-Df_raw[Indexes_train]
      Df_test<-Df_raw[-Indexes_train]
    }else{
      Df_train<-Df_raw[Indexes_train,]
      Df_test<-Df_raw[-Indexes_train,]
    }
    
    l_Epsi_Data_varwise_train[[nameV]]<-Df_train
    l_Epsi_Data_varwise_test[[nameV]]<-Df_test
    
    ### Output to predict: previous time series
    #######
    output_raw<-l_output[[nameV]]
    L_end<-nrow(output_raw)-1
    ### Truncation--> take the observations corresponding to the input residuals
    if(is.vector(Df_raw)){
      output_raw<-output_raw[c(1:L_end)]
      l_output_train[[nameV]]<-output_raw[Indexes_train]
      l_output_test[[nameV]]<-output_raw[-Indexes_train]
    }else{
      output_raw<-output_raw[c(1:L_end),]
      l_output_train[[nameV]]<-output_raw[Indexes_train,]
      l_output_test[[nameV]]<-output_raw[-Indexes_train,]
    }
    
  }
  Ntest<-length(indexes_Kfold[-Indexes_train])
  # Index_delta<-parallel::parSapply(cl = hearts,X =c(1:Ntest),
  #                     FUN = KNN_distancef,f_distance=f_distance,
  #                     K=K_chosen,X_data=l_Epsi_Data_varwise_train,
  #                     X_sim=l_Epsi_Data_varwise_test,hearts=hearts)
  #
  List_Delta<-parallel::parLapply(cl = hearts,X =c(1:Ntest),
                 fun = KNN_distancef_returnMean,
                 f_distance=f_distance,
                 K=K_chosen,X_data=l_Epsi_Data_varwise_train,
                 X_sim=l_Epsi_Data_varwise_test,hearts=hearts,
                Y_data=l_output_train)
  List_Delta<-do.call(c,List_Delta)
  List_Delta_varwise<-list()
  for(Nj in names(l_Epsi_Data_varwise)){
    Inds_founds<-which(names(List_Delta)==Nj)
    Df_delta<-as.data.frame(
      t(rbind(sapply(Inds_founds,FUN = function(x){
      return(List_Delta[[x]])
    }))))
    colnames(Df_delta)<-substr(colnames(Df_delta),
                    start = 2,stop = 3)
    List_Delta_varwise[[Nj]]<-Df_delta
  }
  ### take the selected index and compare the real previous time series (output)
  ### with the predicted one
  Mat_dce<-list()
  for(name_v in names(l_Epsi_Data_varwise)){
    ### Distance per forcing condition
    Data_Xo_v<-l_output_train[[name_v]]
    ### predicted if index
    ### Data_Xo_chosen<-Data_Xo_v[Index_delta,]
    ### predicted if df
    Data_Xo_chosen<-List_Delta_varwise[[name_v]]
    ### real one
    Data_Xo_true<-l_output_test[[name_v]]
    ### Take into account the scale of the forcing condition
    ### use the R^2 coeff
    #Mean_true<-colMeans(Data_Xo_true)
    Weight_summarised<-sum((Data_Xo_true)**2)**(-1/2)
    Gap_mat<-Data_Xo_chosen-Data_Xo_true
    Dces_found<-apply(Gap_mat,MARGIN = 1,
                      FUN = f_distance)
    #Mat_dce[[name_v]]<-Dces_found*Weight_summarised
    Ratio_<-Dces_found*Weight_summarised
    Mat_dce[[name_v]]<-1-Ratio_
    # Mat_dce[[name_v]]<-list("pred"=Data_Xo_chosen,
    #                         "target"=Data_Xo_true)
  }
  ### per obs of the test set--> gives the mean distance
  Overall_Distance<-rowMeans(as.data.frame(Mat_dce))
  return(Overall_Distance)
}
# function(hearts,type_data,list_variable,
#          Name_riskF,link_import_residuals,
#          link_rootVAR_model, K_chosen,
#          f_distance,link_export_generations,
#          cols_gg,prefix_link_resid,
#          subfix_link_Sim,Dates_Johanna,
#          subfix_link_Data,Vect_VARtimes){
#' Title
#'
#' @param hearts 
#' @param f_distance 
#' @param vector_k 
#' @param l_Epsi_Data_varwise 
#' @param l_output 
#' @param number_Kfold 
#'
#' @return
#' @export
#'
#' @examples
K_fold_k_param<-function(hearts,f_distance,
                        number_Kfold,type_data,list_variable,
                        Name_riskF,link_import_residuals,vector_K){

    link_g_chosen<-paste0("residuals_MV/",Name_riskF,"/")
    
    # (1) Import X(M-1) -----------------------------------------------------------
    ######
    repertory<-"data_detrend_Winter/"
    l_Input<-list()
    l_output<-list()
    l_outputF_Epsi_Data<-list()
    d<-length(list_variable)
    l_Epsi_Data_varwise<-list()

    # # (2) Import indexes of exceedances -------------------------------------------
    # #########
    for(name_variable in list_variable){
      if(name_variable=="U"){
        link_data<-paste0(repertory,name_variable,
                             "_trunc_detrend.csv")
        Data<-read.csv(file=link_data)
        
      }else{
        link_data<-paste0(repertory,name_variable,"_detrend.csv")
        Data<-read.csv(file=link_data)
      }
      D<-ncol(Data)
      Data<-Data[,c(2:D)]
      colnames(Data)<-substr(x = colnames(Data),
                               start = 2,stop = 3)
      rownames(Data)<-c(1:nrow(Data))
      l_output[[name_variable]]<-Data

      # Import Resid data -------------------------------------------------------
      Epsi_data<-read.csv(file=paste0(link_import_residuals,
                                      name_variable,"_residuals_",
                                      d,".csv"))
      L2_X<-apply(X = Epsi_data,MARGIN = 1,FUN = calcul_norm_L2)
      l_Input[[name_variable]]<-scale(L2_X)
  }
  N<-nrow(l_Input[[1]])
  Vect_ind_block<-c(1:number_Kfold)
  ### Repartition in K folds for all observations
  indexes_Kfold<-sample(x = Vect_ind_block,size = N,
         replace = TRUE)
  ### Run for all possible values
  l_all_cv<-list()
  for(IND in c(1:length(vector_K))){
    K_j<-vector_K[IND]
    All_CV_j<-t(sapply(X = Vect_ind_block,
          K_fold_k_param_j_test,
          hearts=hearts,f_distance=f_distance,
          l_Epsi_Data_varwise=l_Input,
          l_output=l_output,
          indexes_Kfold=indexes_Kfold,K_chosen=K_j))
    l_all_cv[[IND]]<-do.call(c,All_CV_j)
  }
  return(l_all_cv)
}
