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
                       FUN = calcul_norme_L2)
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
Approach_Angle_One_PCA<-function(LIST_all,Name_for_export,
                                 NbScores_Omega,l_variables,
                                 Indices_exts,d,f_transf){
  Length_T<-ncol(LIST_all[[l_variables[1]]]$transf)
  Omega<-matrix(NA,ncol=d*Length_T,
                nrow =length(Indices_exts))
  
  Ind_final<-0
  LIST_Frechet_OBS<-list()
  L<-length(l_variables)-1
  NE<-l_variables[1]
  d<-length(l_variables)
  for(w in c(1:d)){
    NE<-paste0(Name_for_export,"_",l_variables[w])
  }
  
  for(Z in c(1:d)){
    Beg<-1+(Z-1)*Length_T
    End<-as.numeric(Length_T*(Z))
    nom_v<-l_variables[Z]
    DF_nom<-LIST_all[[nom_v]]$transf[Indices_exts,]
    LIST_Frechet_OBS[[nom_v]]<-DF_nom
    colnames(DF_nom)<-c(1:ncol(DF_nom))
    Excedents_l<-apply(X = DF_nom,MARGIN = 1,
                       FUN = calcul_norme_L2)
    FORME_v<-t(t(DF_nom)%*%diag(Excedents_l^(-1)))
    Omega[,c(Beg:End)]<-f_transf(FORME_v)
  }
  mu_Omega<-colMeans(Omega)
  sigma_Omega<-apply(Omega, 2,FUN = sd)
  Name_for_export_Omega<-paste0(Name_for_export,"_shape_Theta.png")
  png(filename = Name_for_export_Omega,width = 1200,
      height = 600)
  matplot(t(Omega[1:50,]),type="l")
  dev.off()
  ANALYSE_PCA<-FactoMineR::PCA(X =Omega ,
                               scale.unit = TRUE,graph = TRUE,
                               ncp = 10)
  val_lambda<-ANALYSE_PCA$eig[,1]
  vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
  vecteur_propvarexp<-1-vect_diff
  
  Name_for_export_inertia<-paste0(Name_for_export,
                                  "Theta/evol_inertia_Single_Theta")
  
  for(name_V in l_variables[c(1:L)]){
    Name_for_export_inertia<-paste0(Name_for_export_inertia,"_",name_V)
  }
  Name_for_export_inertia<-paste0(Name_for_export_inertia,"_",l_variables[length(l_variables)],
                                  ".png")                               
  png(filename = Name_for_export_inertia,width = 1200,
      height = 600)
  plot(vecteur_propvarexp,type="o",
       ylab="Proportion of unexplained inertia",
       xlab="Number of components")
  abline(v=NbScores_Omega+1,col="red")
  dev.off()
  F_propres<-ANALYSE_PCA$svd$V
  Coordonnees<-ANALYSE_PCA$ind$coord
  Name_for_export_eigen_functions<-paste0(NE,
                                          "_eigen_functions.png")
  png(filename = Name_for_export_eigen_functions,width = 1200,
      height = 600)
  END<-37*d
  plot(c(1:END),ANALYSE_PCA$svd$V[,1],ylab=expression(nu),xlab="Time/variable")
  points(c(1:END),ANALYSE_PCA$svd$V[,2],col="green")
  points(c(1:END),ANALYSE_PCA$svd$V[,3],col="red")
  abline(v=37,col="blue")
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
              "Length_TS"=Length_T))
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
  #exp (f_inv) for inverse transformation
  return(f_transf_inv(Shape_Omega_simul))
}

Inner_k_l<-function(l_kl,Mat){
  k<-l_kl[1]
  l<-l_kl[2]
  Vk<-Mat[,k]
  Vl<-Mat[,l]
  Prod_scal<-as.numeric(t(Vk)%*%Vl)/length(Vk)
  return(Prod_scal)
}
#' Extreme_cov_per_K
#'
#' @param j 
#' @param liste_MV_simul: list[df]. List of time
#' series per tidal cycle
#' @param l_name: list[str]. List of variable names.
#' @param L: int. Number of measures per ts
#' @param d: int. Number of variables.
#' @param Thresh_lg: float. Threshold used for gol

#' @return
#' @export
#'
#' @examples
Extreme_cov_per_K<-function(j,liste_MV_simul,l_name,L,d,
                             Thresh_lg){
  
  Mat<-matrix(NA,nrow = L,
              ncol = d)
  for(l in c(1:d)){
    Name<-l_name[l]
    Df<-liste_MV_simul[[Name]]
    Mat[,l]<-unlist(Df[j,])/Thresh_lg
  }
  ###  NL2_per_var
  Squared_norms<-apply(X = Mat,MARGIN = 2,
        FUN = calcul_norme_L2)**2
  
  ### Scalar_product
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
  Scal_product_XY<-sapply(vect_pair,Inner_k_l,
                Mat=Mat)
  return(c(Scal_product_XY,Squared_norms))
}

#' Launch_extrem_cov_per_K
#'
#' @param liste_MV_Orig 
#' @param l_name 
#' @param k 
#' @param vect_lg 
#' 
#' @return
#' @export
#'
#' @examples
Launch_extreme_cov_per_K<-function(liste_MV_Orig,l_name,k,
                                  vect_lg){
  Order_lg<-order(vect_lg,
                      decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Df<-list()
  for(nameV in l_name){
    Sub_Mat<-liste_MV_Orig[[nameV]][sub_order,]
    Df[[nameV]]<-Sub_Mat
  }
  Dims<-dim(Sub_Mat)
  d<-length(l_name)
  L<-Dims[2]
  All_cov<-as.data.frame(t(sapply(c(1:k),
      FUN = Extreme_cov_per_K,
                  liste_MV_simul = Df,
                  l_name = l_name,L = L,d = d,
                  Thresh_lg=Threshold_lg)))
  Cplmt_cols<-sapply(l_name,function(x){
    return(paste0(x,"norm"))
  })
  colnames(All_cov)<-c("Scal_product",
                          Cplmt_cols)
  END_index<-ncol(All_cov)
  ### Retrieve mean_scal_product
  Sig_XY<-mean(All_cov[,1])
  ### Retrieve product of means of squared norms.
  Denominator_rohXY<-prod(colMeans(All_cov[,c(2:END_index)])**(1/2))
  RHO_XY<-Sig_XY/Denominator_rohXY
  return(list("Sig_XY"=Sig_XY,"RHO_XY"=RHO_XY))
}

Extreme_cov_evol<-function(liste_MV_Orig,l_name,vector_k,
                            Ref_RiskF){
  All_results<-t(sapply(X = vector_k,Launch_extreme_cov_per_K,
                      l_name=l_name,
                      liste_MV_Orig=liste_MV_Orig,
                      vect_lg=Ref_RiskF))
  All_results_df<-data.frame(apply(X = All_results,
                   MARGIN = 2,FUN = unlist))
  return(All_results_df)
}
### Computing confidence band for this statistic
Estimator_rho_sig_knowing_K<-function(indexes_used,list_obs_exts,
                                      l_name,L,d,Thresh_lg){
  Stats_sigma<-as.data.frame(t(sapply(indexes_used,
              Extreme_cov_per_K,
              liste_MV_simul=list_obs_exts,
              l_name=l_name,
              L=L,d=d,
              Thresh_lg=Thresh_lg)))
  Cplmt_cols<-sapply(l_name,function(x){
    return(paste0(x,"norm"))
  })
  colnames(Stats_sigma)<-c("Scal_product",
                           Cplmt_cols)
  END_index<-ncol(Stats_sigma)
  ### Retrieve mean_scal_product
  Sig_XY<-mean(Stats_sigma[,1])
  ### Retrieve product of means of squared norms.
  Denominator_rohXY<-prod(
    colMeans(Stats_sigma[,c(2:END_index)])**(1/2))
  RHO_XY<-Sig_XY/Denominator_rohXY
  return(list("Sig_XY"=Sig_XY,"RHO_XY"=RHO_XY))
}

#' Estimator_1boostrap_sample
#'
#' @param list_obs_exts: list[df]. 
#' @param l_name: vector[str]. Names of variables.
#'
#' @return
#' @export
#'
#' @examples
Estimator_1boostrap_sample<-function(list_obs_exts,
                                     l_name,
                                     Thresh_lg){
  
  DIMS<-dim(list_obs_exts[[1]])
  N<-DIMS[1]
  L<-DIMS[2]
  d<-length(l_name)
  new_indexes<-sample(c(1:N),size = N,replace = TRUE)
  Estimator_1_sample<-Estimator_rho_sig_knowing_K(indexes_used = new_indexes,
                              l_name = l_name,L = L,
                              list_obs_exts = list_obs_exts,
                  d = d,Thresh_lg = Thresh_lg)
  return(Estimator_1_sample)
  
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
