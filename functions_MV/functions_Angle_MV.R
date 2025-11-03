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
  for(w in c(1:length(l_variables))){
    NE<-paste0(Name_for_export,"_",l_variables[w])
  }
  for(Z in c(1:length(l_variables))){
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
  matplot(t(Omega[1:100,]),type="l")
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
  abline(v=NbScores_Omega,col="red")
  dev.off()
  F_propres<-ANALYSE_PCA$svd$V
  Coordonnees<-ANALYSE_PCA$ind$coord
  Name_for_export_eigen_functions<-paste0(NE,
                                          "_eigen_functions.png")
  png(filename = Name_for_export_eigen_functions,width = 1200,
      height = 600)
  plot(c(1:74),ANALYSE_PCA$svd$V[,1],ylab=expression(nu),xlab="Time/variable")
  points(c(1:74),ANALYSE_PCA$svd$V[,2],col="green")
  points(c(1:74),ANALYSE_PCA$svd$V[,3],col="red")
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
