### Approach the angle with several PCA basis.
Approach_Angle_Mult_PCA<-function(Indices_exts,root_export,
                                  LIST_all,l_variables,
                                  list_nb_scores,f_transf){
  LIST_shapes<-list()
  Length_T<-ncol(LIST_all[[l_variables[1]]]$transf)
  LIST_Mu<-list()
  LIST_Sig<-list()
  LIST_PCA_functs<-list()
  NbScores_Omega<-0
  MATRICE_SCORES<-matrix(NA,nrow = length(Indices_exts),
                         ncol=10)
  LIST_Frechet_OBS<-list()
  for(Z in c(1:length(l_variables))){
    nom_v<-l_variables[Z]
    DF_nom<-LIST_all[[nom_v]]$transf[Indices_exts,]
    LIST_Frechet_OBS[[nom_v]]<-DF_nom
    colnames(DF_nom)<-c(1:ncol(DF_nom))
    Excedents_l<-apply(X = DF_nom,MARGIN = 1,
                       FUN = calcul_norm_L2)
    FORME_v<-f_transf(t(t(DF_nom)%*%diag(Excedents_l^(-1))))
    ## With log transformation
    Mu_<-colMeans(FORME_v)
    Sig_<-apply(FORME_v,MARGIN = 2,FUN = sd)
    ANALYSE_PCA<-FactoMineR::PCA(X =FORME_v ,
                                 scale.unit = TRUE,graph = TRUE,
                                 ncp = 10)
    val_lambda<-ANALYSE_PCA$eig[,1]
    vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
    vecteur_propvarexp<-1-vect_diff
    
    Name_for_export_inertia<-paste0(root_export,"/Theta/",nom_v,"/evolMultPCA_inertia_Theta_",nom_v,
                                    ".png")
    print(Name_for_export_inertia)
    Nb_scores<-list_nb_scores[[Z]]
    png(filename = Name_for_export_inertia,width = 1200,
        height = 600)
    plot(vecteur_propvarexp,type="o",
         ylab="Proportion of unexplained inertia",
         xlab="Number of components")
    abline(v=Nb_scores+1,col="red")
    dev.off()
    beg<-NbScores_Omega+1
    End<-NbScores_Omega+Nb_scores
    NbScores_Omega<-NbScores_Omega+Nb_scores
    F_propres<-ANALYSE_PCA$svd$V
    Coordonnees<-ANALYSE_PCA$ind$coord
    MATRICE_SCORES[, beg:End]<-Coordonnees[,1:Nb_scores]
    
    LIST_Mu[[nom_v]]<-Mu_
    LIST_Sig[[nom_v]]<-Sig_
    LIST_PCA_functs[[nom_v]]<-F_propres[,1:Nb_scores]
    
  }
  return(list("Scores"=MATRICE_SCORES[,c(1:NbScores_Omega)],
              "Eigen_functions"=LIST_PCA_functs,
              "Mu"=LIST_Mu,
              "Sig"=LIST_Sig,
              "L_Frechet"=LIST_Frechet_OBS,
             "Length_TS"=Length_T,
             "Nb_scores_omega"=NbScores_Omega))
}
### Approach the angle with one PCA basis.
##########################

#' Approach_Angle_One_PCA
#'
#' @param LIST_all list[dataframe]. Transformed time series
#' @param Name_for_export string. Element to complete 
#' the title of the exported file. 
#' @param NbScores_Omega int. Number of PCA
#' eigenvectors used.
#' @param l_variables vector[string]. Variable names.
#' @param Indices_exts vector[int]. Indices of
#' extreme observations.
#' @param d int. Number of variables. 
#' @param f_transf 
#'
#' @return
#' @export
#'
#' @examples
Approach_Angle_One_PCA<-function(LIST_all,Name_for_export,
                                 NbScores_Omega,l_variables,
                                 Indices_exts,d,f_transf){
  
  Length_T<-list()
  Omega<-list()
  Ind_final<-0
  LIST_Frechet_OBS<-list()
  L<-length(l_variables)-1
  NE<-l_variables[1]
  d<-length(l_variables)
  for(w in c(1:d)){
    NE<-paste0(Name_for_export,"_",l_variables[w])
  }
  for(Z in c(1:d)){
    nom_v<-l_variables[Z]
    DF_base<-LIST_all[[nom_v]]$transf
    DF_nom<-DF_base[Indices_exts,]
    LIST_Frechet_OBS[[nom_v]]<-DF_nom
    colnames(DF_nom)<-c(1:ncol(DF_nom))
    Excedents_l<-apply(X = DF_nom,MARGIN = 1,
                       FUN = calcul_norm_L2)
    FORME_v<-t(t(DF_nom)%*%diag(Excedents_l^(-1)))
    Length_T[[nom_v]]<-ncol(FORME_v)
    Omega[[Z]]<-f_transf(FORME_v)
  }
  ### Retrieve matrix form
  Omega<-do.call(cbind.data.frame,
                 Omega)
  print(ncol(Omega))
  mu_Omega<-colMeans(Omega)
  sigma_Omega<-apply(Omega, 2,FUN = sd)
  Name_for_export_Omega<-paste0(Name_for_export,"_shape_Theta.png")
  png(filename = Name_for_export_Omega,width = 1200,
      height = 600)
  matplot(t(Omega[1:50,]),type="l")
  dev.off()
  ANALYSE_PCA<-FactoMineR::PCA(X =Omega ,
                               scale.unit = TRUE,graph = TRUE,
                               ncp = NbScores_Omega)
  val_lambda<-ANALYSE_PCA$eig[,1]
  vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
  vecteur_propvarexp<-1-vect_diff
  Inertia_explained<-vect_diff[NbScores_Omega+1]
  Name_for_export_inertia<-paste0(Name_for_export,
                                  "Theta/evol_inertia_Single_Theta")
  
  for(name_V in l_variables[c(1:L)]){
    Name_for_export_inertia<-paste0(Name_for_export_inertia,"_",name_V)
  }
  Name_for_export_inertia<-paste0(Name_for_export_inertia,"_",l_variables[length(l_variables)],
                                  ".png")                               
  png(filename = Name_for_export_inertia,width = 800,
      height = 600)
  plot(vecteur_propvarexp,type="o",
       ylab="Proportion of unexplained inertia",
       xlab="Number of components",cex.lab = 3)
  abline(v=NbScores_Omega+1,col="red")
  dev.off()
  F_propres<-ANALYSE_PCA$svd$V
  Coordonnees<-ANALYSE_PCA$ind$coord
  Name_for_export_eigen_functions<-paste0(NE,
                                          "_eigen_functions.png")
  png(filename = Name_for_export_eigen_functions,width = 1200,
      height = 600)
  END<-length(Omega)
  plot(c(1:END),ANALYSE_PCA$svd$V[,1],ylab=expression(nu),xlab="Time/variable")
  points(c(1:END),ANALYSE_PCA$svd$V[,2],col="green")
  points(c(1:END),ANALYSE_PCA$svd$V[,3],col="red")
  dev.off()
  ### Select J eigenfunctions.
  Functions_eigen<-F_propres[,1:NbScores_Omega]
  l_fonction<-Functions_eigen
  MATRICE_SCORES<-Coordonnees[,1:NbScores_Omega]
  return(list("Scores"=MATRICE_SCORES,
              "Eigen_functions"=l_fonction,
              "Mu"=mu_Omega,
              "Sig"=sigma_Omega,
              "L_Frechet"=LIST_Frechet_OBS,
              "Length_TS"=Length_T,
              "Omega_target"=Omega,
              "inertia_explained"=Inertia_explained))
}
######### Simulation from one PCA basis.
Simul_Omega_One_PCA_base<-function(list_Mod_One_PCA,
                                   Simul_coords,
                                   NbScores_Omega,
                                   f_transf_inv){
  Shape_Omega_simul<-function_reconstitution_trajectory_std(Vector_coords = Simul_coords,
                                                            Base_functions_p = list_Mod_One_PCA$Eigen_functions,
                                                            NB_dim =NbScores_Omega,
                                                            mu_t =list_Mod_One_PCA$Mu,
                                                            sd_t =list_Mod_One_PCA$Sig)
  return(f_transf_inv(Shape_Omega_simul))
}
######### Simulation with several PCA basis.
Simul_Omega_Mult_PCA_base<-function(list_Mod_Mult_PCA,M,d,list_nb_scores,
                                    Simul_coords,f_transf_inv){
  LIST_Mu<-list_Mod_Mult_PCA[["Mu"]]
  LIST_Sig<-list_Mod_Mult_PCA[["Sig"]]
  LIST_PCA_functs<-list_Mod_Mult_PCA[["Eigen_functions"]]
  Length_T<-list_Mod_Mult_PCA[["Length_TS"]]
  # Simulation of Omega -----------------------------------------------------
  ##############
  Shape_Omega_simul<-matrix(NA,nrow = M,
                            ncol = Length_T*d)
  N_sum<-0
  for(Z in c(1:d)){
    Nb_scores<-list_nb_scores[[Z]]
    beg<- N_sum+1
    End<- N_sum+Nb_scores
    N_sum<-N_sum+Nb_scores
    Beg_t<-1+Length_T*(Z-1)
    END_t<-Length_T*Z
    Shape_Omega_simul[,Beg_t:END_t]<-function_reconstitution_trajectory_std(
      Vector_coords = Simul_coords[,beg:End],
      Base_functions_p = LIST_PCA_functs[[Z]],
      NB_dim =list_nb_scores[[Z]],
      mu_t =LIST_Mu[[Z]],
      sd_t = LIST_Sig[[Z]])
  }
  return(f_transf_inv(Shape_Omega_simul))
}

Inner_k_l<-function(l_kl,Mat){
  k<-l_kl[1]
  l<-l_kl[2]
  First_fcg<-Mat[[k]]
  Second_fcg<-Mat[[l]]
  names_intersected<-intersect(colnames(First_fcg),
                              colnames(Second_fcg))
  Vk<-as.numeric(unlist(First_fcg[,names_intersected]
                        ))
  L2k<-as.numeric(t(Vk)%*%Vk)/length(Vk)
  Vl<-as.numeric(unlist(Second_fcg[,names_intersected]
                        ))
  L2l<-as.numeric(t(Vl)%*%Vl)/length(Vl)
  ## prod scal
  Prod_scal<-as.numeric(t(Vk)%*%Vl)/length(Vk)
  return(c(Prod_scal,L2k,L2l))

}
#' Extreme_cov_per_K
#'
#' @param j 
#' @param liste_MV_simul list[df]. List of time
#' series per tidal cycle
#' @param l_name list[str]. List of variable names.
#' @param L int. Number of measures per ts
#' @param d int. Number of variables.
#' @param Thresh_lg float. Threshold used for gol
#' @param center_GPD boolean or mean value at each time step
#'
#' @return list. 
#' @export
#'
#' @examples
Extreme_cov_per_K<-function(j,liste_MV_simul,l_name,d,
                             Thresh_lg,center_GPD=FALSE){
  Mat<-list()
  L2_found<-list()
  for(l in c(1:d)){
    Name<-l_name[l]
    Df<-liste_MV_simul[[Name]]
    if(!isFALSE(center_GPD)){
      Df<-Df-center_GPD
    }
    Mat[[l]]<-Df[j,]/Thresh_lg
  }
  
  ### Scalar_product
  Combinations<-t(utils::combn(x = c(1:d),
                               m = 2))
  list_IKL<-list()
  for(j in c(1:nrow(Combinations))){
    IKL<-Inner_k_l(l_kl = Combinations[j,],
                   Mat = Mat)
    list_IKL[[j]]<-IKL
  }
  return(list_IKL)
}
Extreme_cov_per_K<-function(j,liste_MV_simul,l_name,d,
                             Thresh_lg,center_GPD=FALSE){
  Mat<-list()
  L2_found<-list()
  for(l in c(1:d)){
    Name<-l_name[l]
    Df<-liste_MV_simul[[Name]]
    if(!isFALSE(center_GPD)){
      Df<-Df-center_GPD
    }
    Mat[[l]]<-Df[j,]/Thresh_lg
  }
  
  ### Scalar_product
  Combinations<-t(utils::combn(x = c(1:d),
                               m = 2))
  list_IKL<-list()
  for(j in c(1:nrow(Combinations))){
    IKL<-Inner_k_l(l_kl = Combinations[j,],
                   Mat = Mat)
    list_IKL[[j]]<-IKL
  }
  return(list_IKL)
}
Scale_norm_per_Ind<-function(j,liste_MV_simul,l_name,
                             center_GPD){
  Mat<-list()
  d<-length(l_name)
  L2_found<-list()
  for(l in c(1:d)){
    Name<-l_name[l]
    Df<-liste_MV_simul[[Name]]
    if(!isFALSE(center_GPD)){
      Df<-Df-center_GPD
    }
    Mat[[l]]<-Df[j,]
  }
  
  ### Scalar_product
  Combinations<-t(utils::combn(x = c(1:d),
                               m = 2))
  list_IKL<-list()
  for(z in c(1:nrow(Combinations))){
    IKL<-Inner_k_l(l_kl = Combinations[z,],
                   Mat = Mat)
    list_IKL[[z]]<-IKL
  }
  return(list_IKL)
}
#' Title
#'
#' @param l_name 
#' @param k 
#' @param vect_lg 
#' @param Nboot 
#' @param alpha_param 
#' @param CPU_hearts 
#' @param list_sig_L2 
#'
#' @return
#' @export
#'
#' @examples
Launch_extreme_confcov_per_K<-function(l_name,k,
                                       vect_lg,
                                       Nboot,alpha_param,
                                       CPU_hearts,list_sig_L2){
  
  Order_lg<-order(vect_lg,
                  decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Resampled_results<-replicate(n = Nboot,
            Result_1sample_ext(list_sig_L2 = list_sig_L2,
                l_name = l_name,k = k,
                Threshold_lg=Threshold_lg,
                Inds_exts=sub_order,
                CPU_hearts=CPU_hearts))
 
  G<-apply(Resampled_results,
           FUN = function(x){
             df<-melt(as.data.frame(x))},MARGIN=2)
  Melt<-do.call(rbind.data.frame,G)
  Ind_sub<-which(Melt$variable=="gamXY")
  result_resamp<-Melt[Ind_sub,"variable"]
  summary_QQ<-as.data.frame(Melt %>% 
    group_by(pair_vars,variable) %>% 
    summarise(Q1=quantile(value,(alpha_param/2)),
              Q2=quantile(value,1-(alpha_param/2)))
    )
  Indexes<-which(summary_QQ$variable!="k")
  summary_QQ<-summary_QQ[Indexes,]
  summary_QQ$k<-rep(k,nrow(summary_QQ))
  return(summary_QQ)
}

#' Title
#'
#' @param list_sig_L2 
#' @param l_name 
#' @param k 
#' @param Threshold_lg 
#' @param Inds_exts 
#' @param CPU_hearts 
#'
#' @return
#' @export
#'
#' @examples
Result_1sample_ext<-function(list_sig_L2,l_name,k,
                     Threshold_lg,Inds_exts,
                     CPU_hearts){
  
  L_ext<-length(Inds_exts)
  Indexes_sampled<-sample(x = c(1:L_ext),
                          size = L_ext,replace = TRUE)
  Pairs_found<-t(utils::combn(x = c(1:length(l_name)),
                              m = 2))
  END_index<-nrow(Pairs_found)
  list_Sig_RHO<-list()
  for(Z in c(1:END_index)){
    Mat_ZSIG<-list_sig_L2[[Z]]
    ## Take extreme individuals
    Sub_matrix_Z<-Mat_ZSIG[Inds_exts,]
    rownames(Sub_matrix_Z)<-c(1:nrow(Sub_matrix_Z))
    ## Resample obs
    Resampled_mat<-Sub_matrix_Z[Indexes_sampled,]
    rownames(Resampled_mat)<-c(1:nrow(Resampled_mat))
    Resampled_mat<-matrix(Resampled_mat,nrow = nrow(Resampled_mat),
                          ncol=ncol(Resampled_mat))
    ###Gam XY computations
    Scal_prod_penalised<-apply(Resampled_mat,MARGIN = 1,
          FUN = function(x){
            L1<-sqrt(x[2])
            L2<-sqrt(x[3])
            Scal_prod<-x[1]
            y<-Scal_prod/(L1*L2)
            return(y)
          })
    #Scal_prod_penalised<-Resampled_mat[,1]/(Resampled_mat[,2]*Resampled_mat[,3])
    Gam_XY<-mean(Scal_prod_penalised)

    ## Sigma/ Rho computations
    Mu_whole<-colMeans(Resampled_mat)
    Sig<-Mu_whole[1]/(Threshold_lg**(2))
    prod_sig<-prod(Mu_whole[-1]**(1/2))
    Pair_z<-l_name[Pairs_found[Z,]]
    Category<-paste0("(",Pair_z[1],",",Pair_z[2],")")
    list_Sig_RHO[[Z]]<-c(k,Sig,Mu_whole[1]/prod_sig,
                         Gam_XY,Category)
  }
  DF_whole_NEW<-do.call(rbind.data.frame,list_Sig_RHO)
  colnames(DF_whole_NEW)<-c("k","sigXY","rhoXY",
                            "gamXY","pair_vars")
  Toconvert<-c("sigXY","k","rhoXY","gamXY")
  DF_whole_NEW[,Toconvert]<-sapply(Toconvert,
                               function(x)
                               {return(as.numeric(DF_whole_NEW[,x]))
                               })
  return(DF_whole_NEW)
}

#' Title
#'
#' @param liste_MV_Orig 
#' @param l_name 
#' @param k 
#' @param vect_lg 
#' @param center 
#'
#' @return
#' @export
#'
#' @examples
Launch_extreme_cov_per_K<-function(liste_MV_Orig,l_name,k,
                                   vect_lg,center){
  
  Order_lg<-order(vect_lg,
                  decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Df<-list()
  for(nameV in l_name){
    Whole<-liste_MV_Orig[[nameV]]
    Sub_Mat<-Whole[sub_order,]
    Df[[nameV]]<-Sub_Mat
  }
  Dims<-dim(Sub_Mat)
  d<-length(l_name)
  L<-Dims[2]
  All_results<-lapply(c(1:k),
                      FUN = Extreme_cov_per_K,
                      liste_MV_simul = Df,
                      l_name = l_name,d = d,
                      Thresh_lg=Threshold_lg,
                      center_GPD=center)
  Pairs_found<-t(utils::combn(x = c(1:length(l_name)),
                              m = 2))
  END_index<-nrow(Pairs_found)
  list_Sig_RHO<-list()
  for(Z in c(1:END_index)){
    Conv_rho_sig_PairZ<-t(sapply(X = All_results,
                                 FUN = function(x){
                                   return(x[[Z]])
                                 }))
    Mu_whole<-colMeans(Conv_rho_sig_PairZ)
    ##Gam 
    Vgam<-apply(Conv_rho_sig_PairZ,MARGIN = 1,
                               FUN = function(x){
                                 L1<-sqrt(x[2])
                                 L2<-sqrt(x[3])
                                 Scal_prod<-x[1]
                                 y<-Scal_prod/(L1*L2)
                                 return(y)
                               })
    ### remove the weight of the threshold
    Gam<-mean(Vgam)
    prod_sig<-prod(Mu_whole[-1]**(1/2))
    Pair_z<-l_name[Pairs_found[Z,]]
    Category<-paste0("(",Pair_z[1],",",Pair_z[2],")")
    Result_Z<-c(k, Mu_whole[1],Mu_whole[1]/prod_sig,
                Gam,Category)
    list_Sig_RHO[[Z]]<-Result_Z
  }
  DF_whole<-do.call(rbind.data.frame,list_Sig_RHO)
  colnames(DF_whole)<-c("k","sigXY","rhoXY",
                        "gamXY","pair_vars")
  Toconvert<-c("sigXY","k","rhoXY","gamXY")
  DF_whole[,Toconvert]<-sapply(Toconvert,
                               function(x)
                               {return(as.numeric(DF_whole[,x]))
                               })
  return(DF_whole)
}

#' Launch_MVconvergence_per_K
#'
#' @param liste_MV_Orig: list[df]. 
#' @param l_name: vect[str].
#' @param k: int. Number of extreme multivariate extreme time series. 
#' @param vect_lg: vect[float]. Value of the risk function (gol). 
#'
#' @return list. Value for each forcing condition 
#' for a given k of the mean absolute coordinate for several basis functions.
#' @export
#'
#' @examples
Launch_MVconvergence_per_K<-function(liste_MV_Orig,l_name,k,
                                   vect_lg){
  
  Order_lg<-order(vect_lg,
                  decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Df<-list()
  list_conv<-list()
  for(nameV in l_name){
    Sub_Mat<-liste_MV_Orig[[nameV]][sub_order,]
    L2_extj<-apply(X = Sub_Mat,MARGIN = 1,
                   FUN = calcul_norm_L2)
    Angle_extj<-t(t(Sub_Mat)%*%diag(L2_extj^(-1)))
    list_conv[[nameV]]<-Function_conv_univ(
      Shape_d =Angle_extj )
  }
  return(list_conv)
}

#' Title
#'
#' @param liste_MV_Orig 
#' @param l_name 
#' @param vector_k 
#' @param Ref_RiskF 
#'
#' @return Determine the evolution of the sigma and
#'  rho coreelation coefficient
#' @export
#'
#' @examples
Extreme_cov_evol<-function(liste_MV_Orig,l_name,vector_k,
                            Ref_RiskF,CPU_hearts,center=FALSE){
  All_results<-lapply(vector_k,Launch_extreme_cov_per_K,
                      l_name=l_name,
                      liste_MV_Orig=liste_MV_Orig,
                      vect_lg=Ref_RiskF,
                      center=center)
  return(do.call(what = rbind.data.frame,
                 All_results))
}
#' Extreme_cov_confevol
#'
#' @param liste_MV_Orig 
#' @param l_name 
#' @param vector_k 
#' @param Ref_RiskF 
#' @param CPU_hearts object from the parallel package
#' @param center Boolean or float. 
#' @param Nboot int. Number of bootstrap samples.
#' @param alpha_param float. Chosen parameter 
#' of the confidence level
#'
#' @return list. Evolution of the bootstrap confidence 
#' bands with the number of exceedances for each pair
#' of variables.
#' @export
#'
#' @examples
Extreme_cov_confevol<-function(liste_MV_Orig,l_name,vector_k,
                           Ref_RiskF,CPU_hearts,center=FALSE,
                           Nboot,alpha_param){
  
  N<-nrow(liste_MV_Orig[[1]])
  Results<-parLapply(cl = CPU_hearts,X = c(1:N),
         fun = Scale_norm_per_Ind,liste_MV_simul = liste_MV_Orig,
                            l_name = l_name,
                            center_GPD = center
                            )
  Pairs_found<-t(utils::combn(x = c(1:length(l_name)),
                              m = 2))
  END_index<-nrow(Pairs_found)
  list_Sig_RHO<-list()
  for(Z in c(1:END_index)){
    Conv_rho_sig_PairZ<-t(sapply(X = Results,
                                 FUN = function(x){
                                   return(x[[Z]])
                                 }))
    list_Sig_RHO[[Z]]<-Conv_rho_sig_PairZ
  }

  All_results<-lapply(X = vector_k,Launch_extreme_confcov_per_K,
                        l_name=l_name,vect_lg=Ref_RiskF,
                        Nboot=Nboot,
                        alpha_param=alpha_param,
                      CPU_hearts=CPU_hearts,list_sig_L2=list_Sig_RHO)
  return(do.call(what = rbind.data.frame,
                 All_results))
}
#' Convgce_Angle_evol
#'
#' @param liste_MV_Orig list[str: dataframe]. Multivariate
#' time series where the key corresponds to the variable name.
#' 
#' @param l_name vector[str]. Name of forcing conditions
#' @param vector_k vector[int]. Possible number 
#' of exceedances.
#' @param Ref_RiskF vector[float]. Values of the 
#' compound risk function.
#'
#' @return Determine the evolution of the convergence
#' of the angular component for each forcing condition. 
#' @export
#'
#' @examples
Convgce_Angle_evol<-function(liste_MV_Orig,l_name,vector_k,
                            Ref_RiskF){
  
  All_results<-lapply(X = vector_k,
                  Launch_MVconvergence_per_K,
                        l_name=l_name,
                        liste_MV_Orig=liste_MV_Orig,
                        vect_lg=Ref_RiskF)

  return(All_results)
}
  
### Computing confidence band for this statistic
#############

#' Estimator_rho_sig_knowing_K
#'
#' @param list_inputs list[df]. Multivariate time series
#' @param Thresh_lg float. Threshold used.
#' @param CPU_hearts object parallel. CPU used for parallel computing.
#' @param center Bool
#' @return Estimator for the given multivariate time series of the extremal correlation
#'  coefficient
#' @export
#'
#' @examples
Estimator_rho_sig_knowing_K<-function(list_inputs,Thresh_lg,
                                      CPU_hearts,center=FALSE){
  N<-nrow(list_inputs[[1]])
  d<-length(list_inputs)
  l_name<-names(list_inputs)
  All_results<-parLapply(cl = CPU_hearts,X=c(1:N),
       fun = Extreme_cov_per_K,
       liste_MV_simul = list_inputs,
       l_name = l_name,d = d,
       Thresh_lg = Thresh_lg,
       center_GPD=center)
  Pairs_found<-t(utils::combn(x = c(1:length(l_name)),
                              m = 2))
  END_index<-nrow(Pairs_found)
  list_Sig_RHO<-list()
  for(Z in c(1:END_index)){
    Conv_rho_sig_PairZ<-t(sapply(X = All_results,
                                 FUN = function(x){
                                   return(x[[Z]])
                                 }))
    Mu_whole<-colMeans(Conv_rho_sig_PairZ)
    prod_sig<-prod(Mu_whole[-1]**(1/2))
    Pair_z<-l_name[Pairs_found[Z,]]
    Category<-paste0("(",Pair_z[1],",",Pair_z[2],")")
    Vect_input<-c(N, Mu_whole[1],Mu_whole[1]/prod_sig,
                  Category)
    list_Sig_RHO[[Z]]<-Vect_input
    
  }
  Estimator_<-do.call(rbind.data.frame,list_Sig_RHO)
  colnames(Estimator_)<-c("k","sigXY","rhoXY","pair_vars")
  Toconvert<-c("sigXY","k","rhoXY")
  Estimator_[,Toconvert]<-sapply(Toconvert,
             function(x)
             {return(as.numeric(Estimator_[,x]))
             })
  return(Estimator_)
}

#' Estimator_1boostrap_sample
#'
#' @param list_obs_exts list[df]. 
#' @param l_name vector[str]. Names of variables.
#' @param Thresh_lg float. Threshold used.
#' @param CPU_hearts object from the parallel 
#' package to run parallel computation
#' @param center 
#'
#' @return Run the function Estimator_rho_sig_knowing_K
#' for resampled observations
#' @export
#'
#' @examples
Estimator_1boostrap_sample<-function(list_obs_exts,
                                     l_name,Thresh_lg,CPU_hearts,
                                     center){
  DIMS<-dim(list_obs_exts[[1]])
  N<-DIMS[1]
  d<-length(l_name)
  new_indexes<-sample(c(1:N),size = N,replace = TRUE)
  list_sampled<-list()
  for(nameV in l_name){
    Df<-list_obs_exts[[nameV]]
    list_sampled[[nameV]]<-Df[new_indexes,]
  }
  Estimation_1sample<-Estimator_rho_sig_knowing_K(
    list_inputs = list_sampled,
        Thresh_lg = Thresh_lg,
        CPU_hearts = CPU_hearts,
        center=center)
  return(Estimation_1sample)
}

