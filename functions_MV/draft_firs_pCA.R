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
# for(t in c(1,19,37)){
#   for(nameV in l_name){
#     print(nameV)
#     St<-Result_MV$orig[[nameV]][,t]
#     Sub<-Result_MV$resume[[nameV]]
#     sig_t<-Sub$params_transfo[[t]]$sigma
#     gam_t<-Sub$params_transfo[[t]]$mu
#     Kappa_t<-Sub$params_transfo[[t]]$nu
#     # F emp to provide x values -----------------------------------------------
#     logemp_fonction<-Ratio_log_model1(gamma = gam_t,sigma_EGPD = sig_t,
#                                       order_quantiles = Order_quantiles,
#                                       data = St)
# 
#     # Boostrap band -----------------------------------------------------------
#     Ratio_param_bootstrap<-replicate(n = 100,Ratio_generations_EGPD(gamma = gam_t, sigma = sig_t,
#                            kappa =Kappa_t,M = length(St),
#                            order_quantiles = Order_quantiles))
#     Bound_boot_inf<-as.numeric(apply(Ratio_param_bootstrap,MARGIN = 1,FUN = function(x){
#       return(as.numeric(quantile(x,0.025)))}))
#     Bound_boot_sup<-as.numeric(apply(Ratio_param_bootstrap,MARGIN = 1,FUN = function(x){
#       return(as.numeric(quantile(x,0.975)))}))
#     Mean_ratio<-mean(logemp_fonction)
#     df_compar_log_found<-cbind.data.frame(Order_quantiles,
#                                           logemp_fonction,
#                                           Bound_boot_inf,
#                                           Bound_boot_sup)
#     colnames(df_compar_log_found)<-c("q_lev",
#                                      "ratio",
#                                      "bound_inf",
#                                      "bound_sup")
#     GG_kt<-ggplot(data=df_compar_log_found,aes(x=q_lev,y=ratio))+
#       geom_point(aes(col="ratio"))+
#       geom_hline(aes(yintercept = Mean_ratio,col="mean_ratio"))+
#       geom_ribbon(aes(ymin=bound_inf,ymax=bound_sup,linetype="confidence_band"),
#                   alpha=0.02,fill="grey",col="darkblue")+
#       scale_color_manual(values=c("ratio"="blue",
#                                   "mean_ratio"="red",
#                                   "confidence_band"="darkblue"))+
#       scale_linetype_manual("Legend",values=c("confidence_band"=2,
#                                               "mean"=5))+
#       ylab("ratio")+xlab("quantile order")+
#       labs(col="Legend")
#     ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",nameV,"/ComparQ_",nameV,"_t=",t,".png"),
#            plot = GG_kt)
#   }
# 
# }

### Bootstrap confidence band
# ######
# NPY<-length(St)/37
# period_years<-c(0.5,1,2,5,10,20,50,80,
#                 100)
# Qyears<-1-(period_years*NPY)^(-1)
# Qdata<-quantile(St,probs = Qyears)
# QEGPQ<-RL_EGPD_model1(gamma = gam_t,sigma_EGPD = sig_t,
#                       order_quantiles = Qyears,Kappa = Kappa_t)
# RLevel_param_bootstrap<-replicate(n = 100,RLevel_generations_EGPD(gamma = gam_t,
#                                                                   sigma = sig_t,
#                                                                   kappa =Kappa_t,M = length(St),
#                                                                   order_quantiles = Qyears))
# Rl_Bound_boot_mean<-rowMeans(RLevel_param_bootstrap)
# Rl_Bound_boot_inf<-as.numeric(apply(RLevel_param_bootstrap,MARGIN = 1,FUN = function(x){
#   return(as.numeric(quantile(x,0.025)))}))
# Rl_Bound_boot_sup<-as.numeric(apply(RLevel_param_bootstrap,MARGIN = 1,FUN = function(x){
#   return(as.numeric(quantile(x,0.975)))}))
# Rlevel_years<-cbind.data.frame(period_years,Qdata,
#                                Rl_Bound_boot_mean,Rl_Bound_boot_inf,
#                                Rl_Bound_boot_sup)
# colnames(Rlevel_years)<-c("rperiod",
#                           "dataQ",
#                           "egpdQ","Boundinf",
#                           "Boundsup")
# list_t[[nameV]]<-list("RL"=Rlevel_years,
#                       "QQ"=dataQ)
# 
# }
# list_for_plot_EGGPD[[t]]<-list_t
# }

### From Pareto to original.
# TAU<-0.95
# Sim_l<-sapply(X = c(1:ncol(Sim_l)),
#               FUN = function(x){
#                 result<-rep(NA,nrow(Sim_l))
#                 col_j<-Sim_l[,x]
#                 data_j<-Vect_l_function[,x]
#                 Qj<-quantile(data_j,TAU)
#                 ### MEV produces vectors with Frechet margins.
#                 Unif_<-exp(-(col_j)^(-1))
#                 Inds_exts<-which(Unif_>TAU)
#                 Convert_todata<-Unif_
#                 Sub_non_exts<-Unif_[-Inds_exts]
#                 #return the corresponding data quantile for non extremes.
#                 Convert_todata[-Inds_exts]<-sapply(X = Sub_non_exts,
#                                                    function(x){return(quantile(data_j,x))})
#                 # And use the Pareto tails for the extreme ones.
#                 V<-(Unif_[Inds_exts]-TAU)/(1-TAU)
#                 Convert_todata[Inds_exts]<-Seuil_lprime*(1-V)^(-1)
#                 return(Convert_todata)
#               })
# Dot product --MRV test
Nobs<-nrow(Vectors_HTAIL)
# Rank transformation. F--> Pareto
################
# if((convert_HTAIL)==TRUE){
#   Vectors_HTAIL<-apply(X =Vectors_HTAIL,MARGIN = 2,FUN = function(x){
#     Denom<-Nobs+1-rank(x)
#     return(Nobs/Denom)
#   })
# }
# if(length(I)==1){
#   Name_i<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[I],".png")
#   GGi<-Graphics_estimators_gamma(series = Vectors_HTAIL[,I],
#                                  vect_k =Vect_k,
#                                  Title_graphic =" ",
#                                  dims_elt_text = dims_elt_text,
#                                  y_lims = LIMS_Y)
#   ggsave(filename = Name_i,plot =GGi,width=8,
#          height=6)
# }
# if(length(J)==1){
#   Name_j<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_",NAME_Vars[J],".png")
#   obj<-Graphics_estimators_gamma(series = Vectors_HTAIL[,J],
#                                  vect_k =Vect_k,
#                                  Title_graphic =" ",
#                                  dims_elt_text = dims_elt_text,
#                                  y_lims = LIMS_Y)
#   ggsave(filename = Name_j,plot =obj,width=8,
#          height=6)
# }
# L_i<-length(NAME_Vars[I])
# if(L_i>1){
#   First_i<-NAME_Vars[I][1]
#   for(j in c(2:L_i)){
#     First_i<-paste0(First_i,"_",NAME_Vars[I][j])
#   }
# }else{
#   First_i<-NAME_Vars[I]
# }
# L_j<-length(NAME_Vars[J])
# if(L_j>1){
#   Second_j<-NAME_Vars[J][1]
#   for(j in c(2:L_j)){
#     Second_j<-paste0(Second_j,"_",NAME_Vars[J][j])
#   }
# }else{
#   Second_j<-NAME_Vars[J]
# }
# Name_max<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_max_",
#                  First_i,"_",Second_j,".png")
# Max_Risk_functionals<-apply(X =  Vectors_HTAIL[,c(I,J)],MARGIN = 1,
#                             FUN = max)
# GG_max<-Graphics_estimators_gamma(series = Max_Risk_functionals,
#                                   vect_k =Vect_k,
#                                   Title_graphic =" ",
#                                   dims_elt_text = dims_elt_text,
#                                   y_lims = LIMS_Y)
# ggsave(filename = Name_max,plot = GG_max,width=8,
#        height=6)
# 
# Name_min<-paste0("graphiques_MV/Evol_gamma/Evol_gamma_min_",
#                  First_i,"_",Second_j,".png")
# Min_Risk_functionals<-apply(X =  Vectors_HTAIL[,c(I,J)],MARGIN = 1,
#                             FUN = min)
# GG_min<-Graphics_estimators_gamma(series = Min_Risk_functionals,
#                                   vect_k =Vect_k,
#                                   Title_graphic =" ",
#                                   dims_elt_text = dims_elt_text,
#                                   y_lims = LIMS_Y)
# ggsave(filename = Name_min,plot = GG_min,width=8,
#        height=6)