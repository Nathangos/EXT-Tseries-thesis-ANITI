# 
# ######### Convergence work
# ### 1) convergence of the maximum
# dates_import<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]
# Max_<-apply(Result_MV$df,MARGIN = 1,FUN = max)
# Outils_POT_graphique(series = Max_,
#                      seuil = quantile(Max_,0.95),
#                      Q1 = 0.50,
#                      Q2 = 0.975,
#                      dates =dates_import,
#                     titre_variable = "Max")
# # Select a convenient window.
# Th_weights<-Outils_Threshr_choix(Q1 = 0.94,
#                      Q2 = 0.99,
#                      NT_ths = 50,
#                      variable = Max_,
#                      N_v = 3
#                      )
# # Use it to restrict the research window in 
# # threshr.
# Th_weights$u_vec[6]
# Th_weights
# plot(Th_weights$pred_perf[,1])
# plot(Th_weights$pred_perf[,2])
# plot(Th_weights$pred_perf[,3])
# 


# Initials<-mev::gp.fit(Vect_l_function[,1], 35)
# plot(Initials)
# Theta_gp<-as.numeric(Initials$estimate)
# # mev::fit.extgp(Vect_l_function[,1],
# #                model = 3,method="mle",init=c(0.9,Theta_gp))


# # first option, approx  ---------------------------------------------------
# # Bivariate maximum ------------------------------------------------------------
# Vect_l_function<-Vect_l_function
# 
# Angular density with Theta ---------------------------------------------
# R_trigopolar<-sum(R_star1^2+R_star2^2)^(1/2)
# Theta_found<-acos(R_star1/R_trigopolar)
# 
# plot(density(Theta_found))
# 
# VECT_K<-c(10:300)
# par(mfrow=c(2,2))
# Result_1<-sapply(VECT_K,FUN = Hillish_stat,Xi = R_min[First_pop],
#                  Eta = Theta[First_pop])
# plot(VECT_K,Result_1,type="l")
# Result_2<-sapply(VECT_K,FUN = Hillish_stat,Xi = R[-First_pop],
#        Eta = Theta[-First_pop])
# plot(VECT_K,Result_2,type="l",col="blue")
# 
# Result_1neg<-sapply(VECT_K,FUN = Hillish_stat,Xi = R[First_pop],
#                  Eta = -Theta[First_pop])
# plot(VECT_K,Result_1,type="l",col="blue")
# Result_2neg<-sapply(VECT_K,FUN = Hillish_stat,Xi = R[-First_pop],
#                  Eta = -Theta[-First_pop])
# plot(VECT_K,Result_2neg,type="l",col="blue")
# par(mfrow=c(1,1))

# Sy<-summary(Modele_coords)
# Simulations_coord<-VineCopula::RVineSim(N=M,
#                                         RVM = Modele_coords)
# xunif <- seq(0.01, 0.99, length = 100)
# yunif <- seq(0.01, 0.99, length = 100)
# gridcop <- expand.grid(x = xunif, y = yunif)
# Z_line<-which(apply(Modele_coords$Matrix==FALSE,
#             MARGIN = 1,FUN = sum)==0)
# Modele_coords$family
# Prop_inertia<-round(val_lambda/sum(val_lambda),2)*100
# expr_iunif<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
#                                                 j = 1,term = "v")
# 
# expr_junif<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
#                                                 j = 2,term = "v")
# END<-J_Omega-1
# l_couples<-list("1"=c(1,4),"2"=c(1,3),
#                 "3"=c(1,2),"4"=c(4,5))
# SUmy<-summary(Modele_coords)
# for(j_cop in c(1:END)){
#   Couples<-l_couples[[j_cop]]
#   Famij<-SUmy$family[j_cop]
#   PAR<-SUmy$par[j_cop]
#   PAR2<-SUmy$par2[j_cop]
#   Dens_iso<-apply(gridcop,MARGIN = 1,
#                   FUN = function(x){
#                     VineCopula::BiCopPDF(u1 =x[1] ,
#                                          u2 =x[2] ,
#                                          family =Famij ,
#                                          par =PAR ,
#                                          par2 =PAR2 )
#                   })
#   Z_cop<-matrix(Dens_iso,nrow = length(xunif), ncol = length(yunif))
#   contour(xunif, yunif, Z_cop,
#           nlevels = 50,
#           xlab=expr_iunif,
#           ylab=expr_junif,
#           col = "black")
#   points(Unif_coord[,Couples[1]],Unif_coord[,Couples[2]],
#          xlab=expr_iunif,
#          ylab=expr_junif,
#          pch=20,
#          col="blue")
# }