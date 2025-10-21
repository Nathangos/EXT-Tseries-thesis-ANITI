# fnct_centrer_par_colonne<-function(x){
#   return(x-mean(x))
# }
# mu_t<-colMeans(FORME_v)
# SD<-apply(FORME_v,
#           MARGIN = 2,FUN =sd )
# observations_centrees<-apply(X = FORME_v,
#                              FUN = fnct_centrer_par_colonne,
#                              MARGIN = 2)
# nlignes<-nrow(observations_centrees)
# l_Mat_moyenne[[nom_v]]<-mu_t
# l_sigma[[nom_v]]<-SD
# bool_standard<-TRUE
# acp_nom<-PCA(observations_centrees,graph=FALSE,
#              scale.unit =TRUE)
# F_propres<-acp_nom$svd$V
# Coordonnees<-acp_nom$ind$coord
# val_lambda<-acp_nom$eig[,1]
# vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
# vecteur_propvarexp<-1-vect_diff
# plot(vecteur_propvarexp,
#      main=paste0("Evolution de la contribution des axes scale=",
#                  bool_standard," T(",nom_v,")"),type="o",
#      ylab="Proportion de variance non expliquee")
# Nb_scores_nomv<-liste_nb_scores[[nom_v]]
# NB_total<-NB_total+Nb_scores_nomv
# l_fonction[[nom_v]]<-F_propres[,1:Nb_scores_nomv]
# MATRICE_SCORES<-c(MATRICE_SCORES,
#                   c(Coordonnees[,1:Nb_scores_nomv]))
# MS<-matrix(MATRICE_SCORES,byrow = FALSE,nrow=nlignes)
# Unif_coord<-VineCopula::pobs(MS)
# Matrice_C <- c()
# for(j in c(1:NB_total)){
#   lim<-NB_total-1
#   vecteur<-rep(0,NB_total)
#   if(j<lim){
#     vecteur[1:j]<-rep(j,j)
#   }
#   if(j==1){
#     vecteur[1]<-NB_total
#   }
#   if(j==lim){
#     vecteur[1:j]<-rep(j,j)
#     vecteur[1]<-1
#   }
#   if(j==NB_total){
#     vecteur[1]<-lim
#     vecteur[2:NB_total]<-rep(1,lim)
#   }
#   Matrice_C<-c(Matrice_C,vecteur)
# }
# 
# Matrice_C <- matrix(Matrice_C, NB_total, NB_total,byrow = TRUE)


# SImul part --------------------------------------------------------------
# Z<-0
# Lim<-length(liste_nb_scores)-1
# LISTE_shapes<-list()
# vect_found<-c()
# vect_l<-c()
# Coords_ech_orig<-matrix(NA,nrow=M,ncol=NB_total)
# for(indicateur_colonne in c(0:Lim)){
#   nom_s<-names(liste_nb_scores)[indicateur_colonne+1]
#   nb_scores_variable<-liste_nb_scores[[nom_s]]
#   Functions_eigen<-l_fonction[[nom_s]]
#   Mu_t<-l_Mat_moyenne[[nom_s]]
#   Sigma_t<-l_sigma[[nom_s]]
#   init_scores<-Z+1
#   end_scores<-Z+nb_scores_variable
#   Coord_Var_j<-Coords_ech_orig[,c(init_scores:end_scores)]
#   Shape_ACP<-function_reconstitution_trajectory_std(Vector_coords = Coord_Var_j,
#                                                     Base_functions_p = Functions_eigen,
#                                                     NB_dim =nb_scores_variable,
#                                                     mu_t =Mu_t,
#                                                     sd_t = Sigma_t)
#   L2_shape<-apply(X = Shape_ACP,MARGIN = 1,
#                   FUN = calcul_norme_L2)
#   Shape_ACP<-t(t(Shape_ACP)%*%diag(L2_shape^(-1)))
#   LISTE_shapes[[nom_s]]<-Shape_ACP
#   Z<-Z+nb_scores_variable
# }
# for (nom_s in liste_variable){
#   Shape_noms<-LISTE_shapes[[nom_s]]
# }
