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
# if(Approach=="Gauss"){
#   ### To do: correct Sample_cond_g
#   Sim_l<-t(replicate(n =N_sim,
#                      expr = Sample_cond_g_Rjection_Sampling(Mu_vector = Mu_vector,
#                                                             Cov_mat = Cov_mat,
#                                                             g=List_Params_RF[["function"]],
#                                                             Th_g=Mu_Gauss)))
#   ### Conversion to uniform margins.
#   print("here")
#   Sim_l_tf<-apply(Sim_l,MARGIN=2,
#                   FUN = pnorm)
#   ### Mixture model
#   
# }
# if(Approach=="Laplace"){
#   Sim_l<-t(replicate(n =N_sim,
#                      expr = Sample_cond_Laplacian_Rjection_Sampling(
#                        Mu_ALD = Model_ALD$center,
#                        Scatter_ALD=Model_ALD$Scatter,
#                        g=List_Params_RF[["function"]],
#                        Th_g=Mu_Laplace)))
#   Sim_l_tf<-apply(Sim_l,MARGIN=2,
#                   FUN = L1pack::plaplace)
# }
# if(Approach=="HGD"){
#   Sim_l<-t(replicate(n =N_sim,
#                      expr = Sample_cond_HGD_Rjection_Sampling(
#                        Model_HGD = Model_HGD,
#                        g=List_Params_RF[["function"]],
#                        Th_g=Mu_Gauss)))
#   Sim_l_tf<-apply(Sim_l,MARGIN=2,
#                   FUN =pnorm)
# }