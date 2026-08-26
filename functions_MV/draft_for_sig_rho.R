# Bootstrap_all_results_simul<-lapply(X=c(1:Nboot),
#                               function(x){
#                                 Estimator_1boostrap_sample(
#                   list_obs_exts = Simul_for_criterion,
#                   l_name = l_name,
#                   Thresh_lg = Threshold_gamma_0,CPU_hearts = coeurs)
#                                 
#                               })
# M_stat_simul<-apply(X =Bootstrap_all_results_simul,MARGIN = 2,
#               FUN = unlist)
# Bounds_simul<-apply(X = M_stat_simul,FUN=F_bounds_per_column,
#                    MARGIN = 2,Alpha_q=Alpha_param)