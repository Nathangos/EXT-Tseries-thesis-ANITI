rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

#### Import des packages. #####

require(dplyr)
require(ggplot2)
require(reshape2)
require(FactoMineR)
require(extRemes)
require(fda)
require(parallel)
require(VineCopula)
require(grid)
require(gridExtra)
require(mvPot)
require(ggside)
require(mev)
require(scales)
require(gamlss)

### Source des fonctions. ####
#Functions from univariate case. 
source("../fonctions/fonctions_Pareto.R")
source("../fonctions/fonctions.r")
source("../fonctions/fonctions_r_pareto_EVOLTAU.r")

#Specific multivariate functions. 
source("functions_MV/function_simul_MV.r")
source("functions_MV/F_simulations_cond.R")
source("functions_MV/GenericEGPDg_carrer.R")
# import G pareto ---------------------------------------------------------
source("functions_MV/function_estim_simul_gPareto.R")

# Consistency results -----------------------------------------------------
source("functions_MV/Fonctions_consistency.R")
F_names<-ls()
nb_coeurs<-detectCores()-4

# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

# Lancement fonction ------------------------------------------------------
# ------------------------------------------------------------------------

l_root<-"residuals_MV/"
l_name<-c("U","Surcote")
TE<-"residus"
STD_donnees<-TRUE

# Import thresholds -------------------------------------------------------
##############
cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue",
         "model"="red","Storm"="black", "exceedances"="red")
SHS<-rep(0.10,37)
Mat_weights<-cbind(rep(0.10,37),SHS,SHS)
LPU<-Create_list_fromALL(l_name = l_name, 
                         obj_source = Mat_weights)
PL2<-0.05
opt_Frech<-TRUE
# l_nb_SCORES<-Create_list_fromALL(l_name = l_name, 
#                                  obj_source = c(3,3))
M<-2000
show_EGPD<-TRUE
file_dates<-"../ss_tend/HIVER/dates_prises.csv"

# Running code ------------------------------------------------------------
#######

## step length is useful to aid convergence
## if the parameter has a fully param model.
## not theoretically justified if one or more 
## smoothing term.
## https://arxiv.org/html/2404.08331v1, update each parameter 
## with a step length.
list_params_EGPD<-list("n.cyc"=150,"nu.step"=0.1,
                       "tau.step"=0.1,"sigma.step"=0.1,
                       "mu.step"=0.1)
Result_MV<-MarTransfo_TS_EXTGP_MV(type_donnees = "hiver",
                                  coeurs = coeurs,lien_racine = l_root,
                                  liste_noms = l_name,opt_Frech = opt_Frech,
                                  type_entree = TE,p_U = LPU, n.dens=2^(14),
                                  file_dates = file_dates,show_EGPD=show_EGPD,
                                  list_params_EGPD=list_params_EGPD)
Inds_dates_pos<-Result_MV$SUB_set
# db<-as.data.frame(Result_MV$orig$Surcote[,Time_problem])
# colnames(db)<-c("x")
# ## Relaunch
# n.cyc<- list_params_EGPD[["n.cyc"]]
# mu.step<- list_params_EGPD[["mu.step"]]
# sigma.step<- list_params_EGPD[["sigma.step"]]
# nu.step<- list_params_EGPD[["nu.step"]]
# tau.step<- list_params_EGPD[["tau.step"]]
# con<-gamlss.control(n.cyc = n.cyc,mu.step =0.01,
#                     sigma.step = sigma.step, nu.step =nu.step,
#                     tau.step = tau.step,autostep=TRUE)
# Run_one_parameter<- gamlss(x~1, 
#                          data=db,family = EGPD1Family(mu.link = "identity"),
#                          control = con,mu.start=gamma_pb,
#                          sigma.start=sig_pb,nu.start=nu_pb,
#                          method=RS(),sigma.fix = TRUE,
#                          nu.fix = TRUE)
# fitted(Run_one_parameter,"mu",
#        type="response")[1]
Inds_dates_pos<-Result_MV$SUB_set
Result_MV$resume$Surcote$INIT
Sub_set_U<-Result_MV$orig_all$U[Inds_dates_pos,]

# Gof test doable if rotations are forbidden (tawn not work with gof) ------------------------------
#############
First_var<-l_name[1]
Second_var<-l_name[2]
Orig<-Result_MV$orig
ind_var1<-which(names(Orig)==First_var)
ind_var2<-which(names(Orig)==Second_var)

Params1<-Result_MV$resume[[ind_var1]]$params_transfo
print(Params1)
Params2<-Result_MV$resume[[ind_var2]]$params_transfo
print(Params2)
Bool_rotations<-FALSE

# Result of EGPD fitting --------------------------------------------------
################

S_analyse<-Result_MV$orig[[ind_var1]]
hs_analyse<-Result_MV$orig[[ind_var2]]
Order_quantiles<-seq.int(from = 0.1,to = 0.95,length.out = 20)
EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")
Object<-Result_MV$resume$U$Fitting[[1]]
gamlss::rvcov(Object)
muFit <- fitted(Object,"mu")[1]
sigmaFit <- predict(Object,what="sigma", 
                    type="response")[[1]]
nuFit <- predict(Object,what="nu", 
                 type="response")[[1]]
get.K(Object)
# Standard error of each beta
# the conventinal se is too precise 
vcov(Object, type="se")
# the sandwich se is wider  
rvcov(Object, type="se") 
Sigma_bnonrobust<-vcov(Object)
Sigma_b<-rvcov(Object)

summary(Result_MV$resume$Surcote$Fitting[[13]])


# NB_years<-37
# NPY<-length(Inds_dates_pos)/NB_years
# P_lim<-3*NB_years
# periods_years<-c(2,5,10,20,50,80,
#                               100,120,200,400)
# periods_years<-periods_years[which(periods_years<P_lim)]
# P<-1-(NPY*periods_years)^(-1)
# library(stringr)
# list_for_plot_EGGPD<-list()
# Times_chosen<-c(1,13,19,25,31,37)
# list_diag<-list()
# list_egpd<-list()
# list_unit<-list("U"="m/s","Surcote"="m")
# list_ylim_ratio<-list("U"=c(0.1,0.5),
#                       "Surcote"=c(0.5,1.5))
# list_params_EGPD_2<-list("n.cyc"=150,"nu.step"=0.1,
#                        "tau.step"=0.1,"sigma.step"=0.1,
#                        "mu.step"=0.1)
# for(nameV in l_name){
#   MATRIX_EGPD_diag<-matrix(NA,nrow = length(Times_chosen),
#                            ncol=5)
#   Z_line<-1
#   list_t<-list()
#   DATA_EGPD<-matrix(NA,nrow = length(Times_chosen),
#                     ncol=3)
#   for(t in Times_chosen){
#     Fittingname_v<-Result_MV$resume[[nameV]]$Fitting
#     output<-capture.output({
#       plot(Fittingname_v[[t]])})
#     numbers <- str_extract_all(output, "\\d+")
#     numbers2 <- str_extract_all(output, "e-\\d+")
#     numbers3 <- str_extract_all(output, "-")
#     L<-sapply(numbers,function(x){return(length(x))})
#     Sub<-which(L!=0)
#     L2<-sapply(numbers2,function(x){return(length(x))})
#     Sub2<-which(L2!=0)
#     L3<-sapply(numbers3,function(x){return(length(x))})
#     Sub3<-which(L3!=0)
#     Values_found<-sapply(Sub,function(x){
#       elt<-numbers[[x]]
#       return(paste0(elt[1],".",elt[2]))
#     })
#     Values_found2<-rep(" ",length(Sub))
#     real_index<-Sub2-2
#     Values_found2[real_index]<-sapply(Sub2,function(x){
#       power<-numbers2[[x]]
#       return(power)
#     })
#     Sign_value<-rep(" ",length(Sub))
#     real_index<-Sub3-2
#     Sign_value[real_index]<-sapply(Sub3,function(x){
#       sign_v<-numbers3[[x]]
#       return(sign_v)
#     })
#     Values_found<-as.numeric(paste0(Sign_value,
#                                     Values_found,Values_found2))
#     MATRIX_EGPD_diag[Z_line,]<-Values_found
#     St<-Result_MV$orig[[nameV]][,t]
#     Sub<-Result_MV$resume[[nameV]]
#     sig_t<-Sub$params_transfo[[t]]$sigma
#     gam_t<-Sub$params_transfo[[t]]$mu
#     Kappa_t<-(Sub$params_transfo[[t]]$nu)
#     Sigma_parameters<-rvcov(Sub$Fitting[[t]])
#     
#     DATA_EGPD[Z_line,]<-c(Kappa_t,sig_t,
#                           gam_t)
#     Z_line<-Z_line+1
#     Msim<-length(St)
#     Simul_extgp<-mev::rextgp(n=Msim,kappa = (Kappa_t),
#                 sigma = sig_t,xi = gam_t,
#                 type = 1)
#     ### KS confidence band
#     ########
#     Qobject<-qqplot(x = Simul_extgp,y =St )
#     dataQ<-Qobject$qdata
#     GG_plot<-ggplot(dataQ,aes(x=x,y=x))+
#       geom_line(aes(col="model"))+
#       geom_point(aes(y=y,col="data"))+
#       geom_ribbon(aes(ymin=lower,col="confidence_band",
#                       ymax=upper),alpha=0.5,fill="grey",
#                   linetype="dashed")+
#       xlab("Model quantiles")+
#       ylab("Data quantiles")+
#       scale_color_manual(values = cols_)+
#       labs(col="Legend")
#     ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",nameV,"/Q_level_",t,".png"),
#            width=8,height=6)
#     ### Bootstrap confidence band
#     Niv_lower<-1-P
#     Niv_upper<-P
#     ALL_Q<- c(Niv_lower,Niv_upper)
#     
#     N_rep<-100
#     Result_N_rep<-sapply(X =c(1:N_rep), FUN=One_sim,
#                             ALL_Quantiles = ALL_Q,
#                             Kappa_t = Kappa_t,gam_t =gam_t ,
#                             sig_t = sig_t,
#                             list_params_EGPD = list_params_EGPD_2,
#                             Msim = N_rep)
#     Q_025_boot_param<-apply(X = Result_N_rep,MARGIN = 1,
#                             FUN = function(x){return(quantile(x,0.025))})
#     Q_975_boot_param<-apply(X = Result_N_rep,MARGIN = 1,
#                             FUN = function(x){return(quantile(x,0.975))})
#     
#     ## Return level-- upper tail
#     RLevel_em<-sapply(Niv_upper,FUN = function(x){
#       return(as.numeric(quantile(St,probs = x)))
#     })
#     
#     RLevel_egpd<-sapply(Niv_upper,FUN = function(x){
#       return(mev::qextgp(kappa = Kappa_t,
#                          sigma =  sig_t,
#                          xi = gam_t,type = 1,
#                          p=x))
#     })
#     RL_sigma_sq_model<-sapply(Niv_upper,FUN=function(x){
#       Vect_grad<-fprime_Rl_extgp(gam_t = gam_t,
#                              sig_t = sig_t,Kappa_t = Kappa_t,
#                              x)
#       Sigma_sq_rl<-as.numeric(Vect_grad%*%Sigma_parameters%*%Vect_grad)
#       return(Sigma_sq_rl)
#     })
#     print(paste0("t=",t,", variable=",nameV))
#     Rl_sigma_model<-sqrt(RL_sigma_sq_model)
#     Q_sup<-sapply(X = RL_sigma_sq_model,FUN = function(sig){
#       return(qnorm(p = 0.975,mean = 0,sd = sig))
#     })
#     # Bound_sup<-RLevel_egpd+Q_sup
#     # Bound_inf<-pmax(RLevel_egpd-Q_sup,rep(0,length(RLevel_egpd)))
#     Deb<-length(P)+1
#     End<-2*length(P)
#     Bound_inf<-Q_025_boot_param[c(Deb:End)]
#     Bound_sup<-Q_975_boot_param[c(Deb:End)]
#     Bounds<-cbind(Bound_inf,Bound_sup)
#     RL_all<-cbind.data.frame(periods_years,RLevel_em,RLevel_egpd,
#                              Bounds)
#     colnames(RL_all)<-c("period","quantile_emp","quantile_egpd",
#                         "bound_inf","bound_sup")
#     GG_RL<-ggplot(data=RL_all,aes(x=period,y=quantile_emp,col="data"))+
#       geom_point()+
#       geom_line()+
#       geom_point(data=RL_all,aes(x=period,y=quantile_egpd,col="model"),pch=17)+
#       geom_line(data=RL_all,aes(x=period,y=quantile_egpd,col="model"))+
#       scale_x_continuous(labels = custom_labels,transform = "log")+
#       scale_y_continuous(transform="log",
#                        labels = label_number(accuracy = 0.001))+
#       geom_ribbon(mapping = aes(ymin=bound_inf,ymax=bound_sup,col="confidence_band"),alpha=0.15,
#                   fill="grey", linetype = "dashed")+
#       xlab("period (years)")+
#       ylab(paste0("value (",list_unit[[nameV]],")"))+
#       scale_color_manual(values = cols_)+
#       labs(col="Legend")
#     ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",nameV,"/RL_gg_",t,".png"),
#            plot = GG_RL,width=8,height=6)
#     
#     ## Return level-- lower tail
#     RLevel_em_lower<-sapply(Niv_lower,FUN = function(x){
#       return(as.numeric(quantile(St,probs = x)))
#     })
#     
#     RLevel_egpd_lower<-sapply(Niv_lower,FUN = function(x){
#       return(mev::qextgp(kappa = Kappa_t,
#                          sigma =  sig_t,
#                          xi = gam_t,type = 1,
#                          p=x))
#     })
#     pow_behavior<-sapply(Niv_lower,FUN = function(x){
#       return(log(x))
#     })
#     RATIO_pow<-log(RLevel_em_lower)/pow_behavior
#     Bound_inf_lower<-Q_025_boot_param[c(1:length(P))]
#     Bound_sup_lower<-Q_975_boot_param[c(1:length(P))]
#     RL_all_lower<-cbind.data.frame(periods_years,RLevel_em_lower,
#                RLevel_egpd_lower,Bound_inf_lower,Bound_sup_lower)
#     colnames(RL_all_lower)<-c("period","quantile_emp",
#                       "quantile_egpd","bound_inf","bound_sup")
#     DF_log_emp<-cbind.data.frame(periods_years, RATIO_pow)
#     colnames(DF_log_emp)<-c("period","ratio_pow")
#     GG_log_ratio<-ggplot(data=DF_log_emp,aes(x=period,y=ratio_pow,
#                                              col="pow"))+
#       geom_line()+
#       geom_point()+
#       scale_y_continuous(limits = list_ylim_ratio[[nameV]])
#     ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",nameV,"/ratio_log_",t,"_lower.png"),
#            plot = GG_log_ratio,width=8,height=6)
#     GG_RL_lower<-ggplot(data=RL_all_lower,aes(x=period,y=quantile_emp,col="data"))+
#       geom_point()+
#       geom_line()+
#       geom_point(data=RL_all_lower,aes(x=period,y=quantile_egpd,
#                           col="model"),pch=17)+
#       geom_line(data=RL_all_lower,aes(x=period,y=quantile_egpd,
#                           col="model"))+
#       scale_y_continuous(transform="log",
#                          labels = label_number(accuracy = 0.001))+
#       geom_ribbon(mapping = aes(ymin=bound_inf,ymax=bound_sup,col="confidence_band"),
#                   alpha=0.15,fill="grey", linetype = "dashed")+
#       xlab("period (years)")+
#       ylab(paste0("value (",list_unit[[nameV]],")"))+
#       scale_color_manual(values = cols_)+
#       labs(col="Legend")
#     ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",nameV,"/RL_gg_",t,"_lower.png"),
#            plot = GG_RL_lower,width=8,height=6)
#   }
#   DATA_EGPD<-as.data.frame(DATA_EGPD)
#   colnames(DATA_EGPD)<-c("Kappa","sigma","gamma")
#   list_egpd[[nameV]]<-DATA_EGPD
#   MATRIX_EGPD_diag<-as.data.frame(MATRIX_EGPD_diag)
#   colnames(MATRIX_EGPD_diag)<-c("mean","variance",
#                                 "skew_coeff",
#                                 "kurtosis_coeff",
#                                 "Filliben corr coeff")
#   list_diag[[nameV]]<-MATRIX_EGPD_diag
# }
# Fittingname_v[[37]]
# M1<-melt(t(list_diag$U))
# M2<-melt(t(list_diag$Surcote))
# Melting<-rbind.data.frame(M1,M2)
# colnames(Melting)<-c("statistic","time","value")
# Melting$time<-Times_chosen[Melting$time]
# Melting$variable<-c(rep(l_name[1],nrow(M1)),
#                     rep(l_name[2],nrow(M2)))
# Melting$target_normal<-ifelse(Melting$statistic=="mean",
#                               0,
#       no =ifelse(Melting$statistic=="variance",1,
#       no=ifelse(Melting$statistic=="skew_coeff",0,
#       no=ifelse(Melting$statistic=="kurtosis_coeff",3,
#       no=ifelse(Melting$statistic=="Filiben_corr_coeff",1,
#                 Inf)))) )
# GG_normal_diags<-ggplot(data=Melting,aes(x=time,y=value,col=variable,
#                                          group=variable))+
#   geom_point()+facet_wrap(~statistic,scales="free_y")+
#   geom_line(aes(x=time,y=target_normal,col="target"))+
#   scale_color_manual(values = c("U"="red","Surcote"="blue",
#                                 "target"="black"))
# ggsave(filename = paste0("graphiques_MV/EGPD_fitting/Normal_diags.png"),
#        width=8,height=6,plot =GG_normal_diags)
# 
# M1<-melt(t(list_egpd$U))
# M2<-melt(t(list_egpd$Surcote))
# Melting<-rbind.data.frame(M1,M2)
# colnames(Melting)<-c("param","time","value")
# Melting$time<-Times_chosen[Melting$time]
# Melting$variable<-c(rep(l_name[1],nrow(M1)),
#                     rep(l_name[2],nrow(M2)))
# GG_egpd_params<-ggplot(data=Melting,aes(x=time,y=value,col=variable,
#                                          group=variable))+
#   geom_point()+facet_wrap(~param,scales="free_y")+
#   scale_color_manual(values = c("U"="red","Surcote"="blue",
#                                 "target"="black"))
# ggsave(filename = paste0("graphiques_MV/EGPD_fitting/Param_found.png"),
#        width=8,height=6,plot =GG_egpd_params)
#     ### Bootstrap confidence band
#     ######
#     NPY<-length(St)/37
#     period_years<-c(0.5,1,2,5,10,20,50,80,
#                     100)
#     Qyears<-1-(period_years*NPY)^(-1)
#     Qdata<-quantile(St,probs = Qyears)
#     QEGPQ<-RL_EGPD_model1(gamma = gam_t,sigma_EGPD = sig_t,
#                    order_quantiles = Qyears,Kappa = Kappa_t)
#     RLevel_param_bootstrap<-replicate(n = 100,RLevel_generations_EGPD(gamma = gam_t, 
#                             sigma = sig_t,
#                             kappa =Kappa_t,M = length(St),
#                             order_quantiles = Qyears))
#     Rl_Bound_boot_mean<-rowMeans(RLevel_param_bootstrap)
#     Rl_Bound_boot_inf<-as.numeric(apply(RLevel_param_bootstrap,MARGIN = 1,FUN = function(x){
#                              return(as.numeric(quantile(x,0.025)))}))
#     Rl_Bound_boot_sup<-as.numeric(apply(RLevel_param_bootstrap,MARGIN = 1,FUN = function(x){
#                              return(as.numeric(quantile(x,0.975)))}))
#     Rlevel_years<-cbind.data.frame(period_years,Qdata,
#                                    Rl_Bound_boot_mean,Rl_Bound_boot_inf,
#                                    Rl_Bound_boot_sup)
#     colnames(Rlevel_years)<-c("rperiod",
#                               "dataQ",
#                               "egpdQ","Boundinf",
#                               "Boundsup")
#     list_t[[nameV]]<-list("RL"=Rlevel_years,
#                           "QQ"=dataQ)
#       
#   }
#   list_for_plot_EGGPD[[t]]<-list_t
# }
# Unit_used<-list("m","m")
# names(Unit_used)<-l_name
# Name_chosen<-l_name[2]
# t_chosen<-37
# Object_Rl<-list_for_plot_EGGPD[[t_chosen]][[Name_chosen]]$RL
# GGRL<-ggplot(data=Object_Rl,aes(x=rperiod,y=dataQ))+
#   geom_point(aes(col="data"))+
#   geom_line(aes(y=egpdQ,col="model"))+
#   scale_color_manual(values = cols_)+
#   xlab("Period P (years)")+
#   ylab(paste0("return level (",Unit_used,")"))+
#   geom_ribbon(aes(ymin=Boundinf,ymax=Boundsup,col="confidence_band"),alpha=0.5,fill="grey",
#               linetype="dashed")+
#   scale_y_continuous(transform="log",
#                      labels = label_number(accuracy = 0.1))+
#   scale_x_continuous(labels = custom_labels,transform = "log")+
#   labs(col="Legend")
# ggsave(filename = paste0("graphiques_MV/RL_EGPD",Name_chosen,"_t=",t_chosen,".png"),
#        plot = GGRL,width=8,height=6)
# dev.off()

# Comparing with Frechet distrib ------------------------------------------
###########
Vect_l_function<-Result_MV$df
Name_orig<-"residuals_MV/RiskF"
for(j in c(1:length(l_name))){
  Name_orig<-paste0(Name_orig,"_",l_name[j])
}
write.csv(x =Vect_l_function,
          file = paste0(Name_orig,".csv"),
          )

# To export, make it easy bro !! ------------------------------------------

Orig<-Result_MV$orig
summary(Vect_l_function)
l_delta<-seq.int(from = 0.1,to = 0.95,
                 length.out = 10)
Frechet_moments<-gamma(1-l_delta)
l_moments<-matrix(NA,nrow = length(l_delta),
                  ncol=ncol(Vect_l_function))
Cst_found<-l_moments
for(z in c(1:length(l_delta))){
  delta<-l_delta[z]
  l_moments[z,]<-colMeans(Vect_l_function^delta)
  Cst_found[z,]<-(l_moments[z,]/Frechet_moments[z])^(1/delta)
}
par(mfrow=c(1,ncol(Vect_l_function)))

# Comparing the evol of the function --------------------------------------
# with the evol for a Frechet distribution --------------------------------

for(nc in c(1:ncol(Vect_l_function))){
  plot(l_delta,l_moments[,nc],main=paste0("Evolution 
     of power moments for the ",nc," th variable"))
  lines(l_delta,Frechet_moments,col="red")
  plot(l_delta,Cst_found[,nc])
}
Cst_random<-Cst_found[5,]
Unif_failed<-exp(-(Vect_l_function/Cst_random)^(-1))
summary(Unif_failed)
Unif_failed<-Unif_failed
for(nc in c(1:ncol(Vect_l_function))){
  print("AD test post Frechet cdf application")
  pval<-goftest::ad.test(x = Unif_failed[,nc],
                         null = "punif")$p.value
  print(pval)
  plot(density(Unif_failed[,nc]),main=" ")
}
# U<--1/log(runif(n=1000))
# expnc<-1/U
z_test<-1
expnc<-(Vect_l_function[,z_test]/Cst_random[z_test])
Theta<-mev::gp.fit(xdat = expnc,threshold = 0)$estimate
Theta
## We do have a positive gamma but lower than 1
All_theta_mle<-mev::fit.extgp(data = expnc,model = 1,
               method="mle",
               init = c(1,Theta),
               )
All_theta_mle$fit

mtext(paste0("Density post Frechet cdf/compare with Unif"),
      outer=TRUE,line=-3)

par(mfrow=c(1,1))
# Evol Poisson aspect -----------------------------------------------------
# dates_import<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[Result_MV$SUB_set,2]
# d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
# Decimal_date<-lubridate::decimal_date(d_POIXCT)
# Nb_years<-37
# x_ <- seq(1, 200, length = 100)
# y_ <- x_
# grid <- expand.grid(x = x_, y = y_)
# 
# Result_Dindex<-apply(X = grid ,FUN = function(x){
#   U<-as.numeric(x[1])
#   V<-as.numeric(x[2])
#   return(DIndex_Poisson_Thresholds(Vect_l_function = Vect_l_function,
#                                    U = U,V = V,
#                                    Decimal_date = Decimal_date,
#                                    Nb_years = 37))},MARGIN=1)
# Result_Dindex
# DIndex<-Result_Dindex
# Little_delta<-as.character(ifelse(abs(DIndex-1)<0.05,1, 0))
# Extreme_aspect<-cbind.data.frame(grid,DIndex,Little_delta)
# colnames(Extreme_aspect)<-c("grid_x","grid_y","DIndex",
#                             "ind_Poisson")
# 
# # ggplot(data=Extreme_aspect,aes(x=grid_x,y=grid_y,z=DIndex))+
# #   geom_tile(aes(fill =DIndex )) +
# #   scale_x_continuous(transform = "log")+
# #   scale_y_continuous(transform = "log")+
# #   stat_contour(aes(fill = ..level..), geom = 'polygon', binwidth = 0.005)+
# #   scale_fill_viridis_c()
# 
# # Results -----------------------------------------------------------------
# 
# require(plotly)
# plot_ly(x=Extreme_aspect$grid_x, 
#         y=Extreme_aspect$grid_y, z=Extreme_aspect$DIndex, type="scatter3d", mode="markers",
#         color = Extreme_aspect$DIndex)


# Asymptotic dependence ? -------------------------------------------------
##########
# For t,t', relation between the variables --------------------------------
##########

Epsi_1<-Orig[[ind_var1]]
Epsi_2<-Orig[[ind_var2]]
Epsi_tf1<-Result_MV$resume[[ind_var1]]$transf
Epsi_tf2<-Result_MV$resume[[ind_var2]]$transf

# Graphical analysis of asymp dependencies --------------------------------
###########
t<-1
t_prime<-37
ROOT_analysis<-paste0(getwd(),'/graphiques_MV/')
ROOT_analysis
Export_chimeas_result<-function(i,j,Orig,t1,t2){
  Epsi_1<-Orig[[i]]
  var1<-names(Orig)[i]
  Epsi_2<-Orig[[j]]
  var2<-names(Orig)[j]
  png(paste0(ROOT_analysis,"asymp_dependencies/ext_dep_",var1,"_t=",t1,"_",var2 ,"_tprime=",t2,".png")
             ,width = 1200, height = 600)
  par(mfrow = c(1, 2),  # 1 row, 2 plots
      cex.lab = 1,    # axis label size
      cex.axis = 1,   # axis tick label size
      cex.main = 1,   # main title size
      cex = 1)
  Two<-expression(bar(chi))
  POT::chimeas(data = cbind(Epsi_1[,t1],Epsi_2[,t2])
               ,which = 1,ask = FALSE,cex.lab=1.5)
  abline(h=0,col="red")
  POT::chimeas(data = cbind(Epsi_1[,t1],Epsi_2[,t2])
               ,which = 2,ask = FALSE,cex.lab=1.5,
               ylabs = rep(NA,2))
  mtext(expression(bar(chi)),cex=1.5,las=2,
        outer=TRUE,line=-20.5)
  par(mfrow=c(1,1))
  dev.off()
}
I<-2
J<-2
s<-1
t<-37
Export_chimeas_result(i = I,j = J,Orig = Orig,
                      t1 =s,t2 =t )

# Same thing for vector of risk functions ---------------------------------
png(paste0(ROOT_analysis,"asymp_dependencies/riskF_Dep_",
           l_name[I],"_",l_name[J],".png")
    ,width = 1200, height = 600)
par(mfrow = c(1, 2),  # 1 row, 2 plots
    cex.lab = 1,    # axis label size
    cex.axis = 1,   # axis tick label size
    cex.main = 1,   # main title size
    cex = 1)
POT::chimeas(data = cbind(Vect_l_function[,I],
                          Vect_l_function[,J])
             ,which = 1,ask = FALSE,
             cex.lab=1.5)
abline(h=0,col="red")
POT::chimeas(data = cbind(Vect_l_function[,I],
                          Vect_l_function[,J])
             ,which = 2,ask = FALSE,
             cex.lab=1.5,
             ylabs = rep(NA,2))
mtext(expression(bar(chi)),cex=1.5,las=2,
outer=TRUE,line=-20.5)
dev.off()
par(mfrow=c(1,1))

Pval_Test<-POT::tailind.test(data =cbind(Epsi_1[,t],Epsi_2[,t_prime]),emp.trans =TRUE,c=-0.1)
p_nullIndep<-Pval_Test$stats[,2]
print(p_nullIndep)
Message<-ifelse(all(as.numeric(p_nullIndep>0.05)),yes = "Asymp independence kept",
                no = "Asymp independence rejected")
print(Message)

Hnull_asymp_dep<-extRemes::taildep.test(x = Epsi_1[,t],y=Epsi_2[,t_prime],
                                        cthresh = -0.40)
p_nulldep<-as.numeric(Hnull_asymp_dep$p.value)

# Exceedance rate should be between 10 and 15 -----------------------------
#################
print(c(p_nulldep,Hnull_asymp_dep,
        Hnull_asymp_dep$parameter[4]))
Message<-ifelse(p_nulldep>0.05,yes = "Asymp dependence kept",
                no = "Asymp dependence rejected for conf. level 0.05")
print(Message)

# Work on hidden regular variations ---------------------------------------
# Sample in decreasing order ----------------------------------------------
R_star1<-rank(x =Epsi_tf1[,t])
R_star2<-rank(x =Epsi_tf1[,t_prime])
Graphics_estimators_gamma(series = Epsi_tf1[,t],
                          vect_k = c(10:500),
                          Title_graphic = paste0(First_var,"_",
                                                 t))
Graphics_estimators_gamma(series = Epsi_tf2[,t_prime],
                          vect_k = c(10:500),
                          Title_graphic = paste0(Second_var,"_",
                                                 t_prime))
Min_<-apply(cbind(Epsi_tf1[,t],Epsi_tf2[,t_prime]),
            MARGIN = 1,FUN = min)
Graphics_estimators_gamma(series =Min_,
                          vect_k = c(10:500),
                          Title_graphic = paste0("Min_(",First_var,
                          "_",t,",",Second_var,"_",t_prime,")"))

# Stats to confirm/or not HRV ---------------------------------------------
###############
R_min<-apply(cbind(R_star1,R_star2),MARGIN = 1,FUN = min)
Frac1<-R_star1/R_star2
Frac2<-Frac1^(-1)
First_pop<-which(Frac1>Frac2)
Theta<-apply(cbind(Frac1,Frac2),MARGIN = 1,FUN = max)
# Measures of dependence --------------------------------------------------
###############
par(mfrow=c(1,2))
POT::chimeas(data =Result_MV$df,
             which = 1,ask = FALSE)
abline(h=0,col="red")
POT::chimeas(data = Result_MV$df
             ,which = 2,ask = FALSE)
abline(h=0,col="red")
par(mfrow=c(1,1))
for(j in c(1:length(l_name))){
  obj1<-Graphics_estimators_gamma(series = Result_MV$df[,j],
                                  vect_k =c(10:500),
                                  Title_graphic =" ")
  ggsave(filename = paste0("graphiques_MV/EGPD_fitting/",l_name[j],
                           "/shift_tail_",l_name[j],".png")
         ,plot = obj1)
}
# (Hnull) asymptotic independence 
###################

Pval_Test<-POT::tailind.test(data =Result_MV$df,emp.trans =TRUE,c=-0.1)
p_nullIndep<-Pval_Test$stats[,2]
print(p_nullIndep)
Message<-ifelse(all(as.numeric(p_nullIndep>0.05)),yes = "Asymp independence kept",
                no = "Asymp independence rejected")
print(Message)

# (Hnull) Asymptotic dependence
#############
Hnull_asymp_dep<-extRemes::taildep.test(x = Result_MV$df[,1],y=Result_MV$df[,2],
                                        cthresh = -0.35)
p_nulldep<-as.numeric(Hnull_asymp_dep$p.value)

# Exceedance rate should be between 10 and 15 -----------------------------
###############
print(c(p_nulldep,Hnull_asymp_dep,
        Hnull_asymp_dep$parameter[4]))
Message<-ifelse(p_nulldep>0.05,yes = "Asymp dependence kept",
                no = "Asymp dependence rejected")
print(Message)

# Work on null hypothesis -------------------------------------------------
##############
for(NAME_for_test in l_name){
  EXTGP<-Result_MV$resume[[NAME_for_test]]$transf
  ad_<-rep(NA,ncol(EXTGP))
  ks_<-rep(NA,ncol(EXTGP))
  cvm_<-rep(NA,ncol(EXTGP))
  for(t in c(1:ncol(EXTGP))){
    ad_[t]<-goftest::ad.test(x = EXTGP[,t],null = extRemes::"pevd",
                             scale=1,shape=1,loc=1,type="GEV", 
                             threshold=1)$p.value
    ks_[t]<-ks.test(EXTGP[,t],extRemes::"pevd",
                    scale=1,shape=1,loc=1,type="GEV",
                    threshold=1)$p.value
    cvm_[t]<-goftest::cvm.test(EXTGP[,t],extRemes::"pevd",
                               scale=1,shape=1,loc=1,type="GEV",
                               threshold=1)$p.value
  }
  plot(c(1:37),ad_,ylim = c(0,1),ylab="p value gof tests")
  points(c(1:37),ks_,col="red")
  points(c(1:37),cvm_,col="green")
  abline(h=0.05,col="black")
}
Vect_l_function
Vect_l_function<-as.data.frame(Vect_l_function)

# Maximums used -----------------------------------------------------------
# NOT affected by assumption of frechet for indual ------------------------

NPY_per_block_OPTIM<-40
NRows<-nrow(Vect_l_function)
M_max<-round(NRows/NPY_per_block_OPTIM)
BLOCKS<-as.character(round(c(1:NRows)/NPY_per_block_OPTIM))
Calul_max<-cbind.data.frame(Vect_l_function,BLOCKS)
colnames(Calul_max)<-c(l_name,"block")

Vect_input<-Calul_max %>% group_by(block) %>%
summarise(across(l_name,~max(.x,na.rm = TRUE)))

#normalisation
Vect_GEV<-as.data.frame(Vect_input[,1+c(1:length(l_name))])/NPY_per_block_OPTIM

delta_moment<-0.5
value_found<-colMeans(Vect_GEV^delta_moment)
value_std_frechet<-factorial(-delta_moment)
scale_frechet_d<-(value_found/value_std_frechet)^(1/delta_moment)
Result_std<-matrix(NA,nrow = nrow(Vect_GEV),
                   ncol=ncol(Vect_GEV))
Fill<-matrix(NA,ncol=ncol(Vect_GEV),nrow=3)
Result_std<-Vect_GEV
for(j in c(1:ncol(Vect_GEV))){
  Result_std[,j]<-Vect_GEV[,j]/scale_frechet_d[j]
  pval<-goftest::ad.test(Vect_l_function[,j],null = extRemes::"pevd",
                         scale=scale_frechet_d[j],shape=1,
                         loc=scale_frechet_d[j],type="GEV")$p.value
  pval2<-ks.test(Vect_l_function[,j],extRemes::"pevd",
                 scale=scale_frechet_d[j],shape=1,
                 loc=scale_frechet_d[j],type="GEV")$p.value
  pval3<-goftest::cvm.test(Vect_l_function[,j],extRemes::"pevd",
                           scale=scale_frechet_d[j],shape=1,
                           loc=scale_frechet_d[j],type="GEV")$p.value
  Fill[,j]<-c(pval,pval2,pval3)
}
Fill
summary(Vect_GEV)
summary(Vect_l_function)
Result_std<-as.matrix(Result_std)
delta<-0.05
w<-seq.int(0, 1,by=delta)
Anonparam<-evd::abvnonpar(w,data =Result_std)
plot(w,Anonparam,type="l")
TML<-"param"

## Bivariate extreme value model -------------------------------------------
# First option to choose model --------------------------------------------

typemod<-"log"
SYMbool<-FALSE

Mod_biv<-evd::fbvevd(x = Result_std,
                     method="BFGS",model=typemod,
                     sym=SYMbool,std.err = FALSE)
Object<-plot(Mod_biv,which = 4)
Mod_biv$estimate
if(typemod=="log"){
  ALPHA_Mstable<-NULL
}
if(typemod=="bilog"){
  ALPHA_Mstable<-c(Mod_biv$estimate[["alpha"]],
                 Mod_biv$estimate[["beta"]])
}
if(typemod%in%c("alog","hr")){
  ALPHA_Mstable<-Mod_biv$estimate[["dep"]]
}

Liste_transf<-list()
for(name in l_name){
  obj<-Result_MV$resume[[name]]
  Liste_transf[[name]]<-obj$transf
}

# 3) simulation procedure -------------------------------------------------
### ----------------------------------------------------------------------

# Dombry/ribatet wrongly interpretated ------------------------------------

# so be careful with frechet margins assumptions --------------------------
# Better solution with EGPD interface for none_param transformation to exp. -------------------------------------
# EGPD) 0.75, normally 0.95

POT::chimeas(data = Result_MV$df,ask = FALSE)

par(mfrow=c(1,1))

# 1)  ---------------------------------------------------------------------
## Convergence of angle measure of the g Pareto process 
#Name_riskF<-"sum_penalized"
q<-(20)
W_chosen<-c(0.3,0.7)
Name_riskF<-"Weighted_NQ"
Name_riskF2<-Name_riskF
if(Name_riskF=="sum"){
  Risk_f<-function(x){return(sum(x))}
  wFun<-weightFun_sum
  dwFun<-DweightFun_sum
}
if(Name_riskF=="max"){
  Risk_f<-function(x){return(max(x))}
  wFun<-NA
  dwFun<-NA
}
if(Name_riskF=="sum_penalized"){
   Risk_f<-function(x){return(sum_penalized(x = x,vect_w = c(0.7,0.3)))}
   wFun<-function(x,u){
     weightFun_sum_penalized(x = x,u = u,
                             vect_w = W_chosen)
   }
   dwFun<-function(x,u){
     DweightFun_sum_penalized(x = x,u = u,
                              vect_w = W_chosen)
   }
}

if(Name_riskF=="Weighted_NQ"){
  Name_riskF2<-paste0(Name_riskF2,"_q=",q)
  Risk_f<-function(x){
    return(Norm_q(x=x,q=q,
                  weights_nq = c(0.7,0.3)))
  }
  wFun<-function(x,u){
    return(weightFun_nq(x = x,u = u,
                        q = q,
                        weights_nq =c(0.7,0.3)))
  }
  dwFun<-function(x,u){
    return(DweightFun_nq(x = x,u = u,
                         q = q,weights_nq = c(0.7,0.3)))
  }
}
path <- paste0(getwd(),"/residuals_MV")
path_graph<-paste0(getwd(),"/graphiques_MV/resid")
newfolder <- Name_riskF2
try({
  dir.create(file.path(path, newfolder))
})
ROOT_export<-paste0(path,"/",Name_riskF2,"/")
try({
  dir.create(file.path(path_graph, newfolder))
})
ROOT_export_graph<-paste0(path_graph,"/",
                          Name_riskF2,"/")
try({
  TEND_path<-paste0(ROOT_export_graph,"Tend")
  dir.create(file.path(TEND_path))
  Lg_path<-paste0(ROOT_export_graph,"RisKfunctions")
  dir.create(file.path(Lg_path))
  Theta_path<-paste0(ROOT_export_graph,"Theta")
  dir.create(file.path(Theta_path))
  PCA_path<-paste0(ROOT_export_graph,"PCA")
  dir.create(file.path(PCA_path))
  Extremo_path<-paste0(ROOT_export_graph,"extremes")
  dir.create(file.path(Extremo_path))
  
  
})

# 2) Choice of the threshold -------------------------------------------------
Riskf_data<-apply(X = Vect_l_function,MARGIN = 1,FUN = Risk_f)

# Distance depends on a truncation of the density -------------------------
# So quite important ! ----------------------------------------------------

TS_param<-0.50
Result_k0_metric1<-tea::mindist(data = Riskf_data,method = "ks",ts =TS_param )
kstar1<-Result_k0_metric1$k0
print(paste0("We would take ",kstar1," extreme MV time series with KS"))
Result_k0_metric2<-tea::mindist(data = Riskf_data,method = "mad",ts =TS_param )
kstar2<-Result_k0_metric2$k0
print(paste0("We would take ",kstar2," extreme MV time series with MAD distance"))

# Evolution of the metric -------------------------------------------------
######

k_metric<-mindist_update(data = Riskf_data,ts = TS_param,method = "ks")
kZero_found<-which.min(k_metric$value_metric)
k_cond<-37
if(!is.null(k_cond)){
  if(kZero_found<k_cond){
    end<-nrow(k_metric)
    sub_kmetric<-k_metric[k_cond:end,]
    kZero_found<-which.min(sub_kmetric$value_metric)+k_cond
  }
}
png(filename = paste0(Lg_path,"/graph_metric_lg.png"),
    width = 500,height = 500)
plot(k_metric$Nb_k,k_metric$value_metric,xlab="Number of exceedances",
     ylab="score_metric",type="l")
abline(v=kZero_found,col="red")
kZero_found
dev.off()
summary(Riskf_data)
PLOT_RISKF<-Graphics_estimators_gamma(series = Riskf_data,
                          vect_k =c(10:150),
                          Title_graphic =" ")
## ~75 individuals, first stability period.
ggsave(filename = paste0(ROOT_export_graph,"Shape_lg_.png"),
       plot = PLOT_RISKF,width = 8,height=6)
#kZero_found<-100
Qfound<-round(1-(kZero_found/length(Riskf_data)),2)
#Qfound<-0.80
Analyse_extreme_proj_gfunction(base_RV = Vect_l_function,
                               L = 400,M1 = kZero_found-50,
                              M2 =kZero_found+50,function_g = Risk_f)
# ### 2-bis) convergence of angle and correlation
Max_nb_exts<-500
vectk<-seq.int(from = 50,
               to = Max_nb_exts,
               by = 1)
Con2results<-sapply(vectk,FUN = Conv_scalar_product,
                       l_name_variables = l_name,
                       list_TS = Liste_transf,g_function=Risk_f)
rho_corr<-unlist(Con2results[1,])
plot(vectk,rho_corr,xlab="Number exceedances",
     ylab="mean_beta",type="l")
plot(vectk,unlist(Con2results[2,])
     ,xlab="Number excedents",
     ylab="var_rho",type="l")
Matrix<-sapply(Con2results[3,],unlist)
par(mfrow=c(2,2))
apply(X = Matrix,FUN=
        function(x,vect_k){
          y<-unlist(x)
          plot(vect_k,y,type="l",
               xlab="Number exceedances",
               ylab="First moment")
        },MARGIN=1,vect_k=vectk)
par(mfrow=c(1,1))

# Choice of extreme functional --------------------------------------------
### Behaviour of pairs of U
source("functions_MV/Hidden_RV_functions.R")
T<-1
S<-2
Var1<-l_name[1]
Var2<-l_name[1]
U<-Result_MV$resume[[Var1]]$transf[,T]
V<-Result_MV$resume[[Var2]]$transf[,S]
Name_vars<-c(paste0(Var1,"_",T),
             paste0(Var2,"_",S))
Vect_k<-c(10:400)
LIMS_Y<-c(0,2)
DIMS_plot<-c(15,15,14,14)
Run_diagnostics_Gamma_G(LIMS_Y =c(-1,1),
                        dims_elt_text = DIMS_plot,
                        Vectors_HTAIL = cbind(U,V),
                        I = 1,J = 2,NAME_Vars = Name_vars,
                        Vect_k = Vect_k,q=0.8)
par(mfrow=c(1,1))
Vect<-cbind(U,V)
Min_uv<-apply(X = Vect,MARGIN = 1,FUN = min)
k<-500
min_tau<-quantile(Min_uv,1-(k/length(Min_uv)))
exceed<-which(Min_uv>min_tau)
Vect_extremes<-Vect[exceed,]

Indexes_<-sapply(X = c(1:nrow(Vect_extremes)),
                 FUN = function(x){
                   return(ifelse(Vect_extremes[x,2]>Vect_extremes[x,1],
                                 yes = 1,no = 0))})
Subset_1<-Vect[which(Indexes_==1),]
Theta_1<-Vect_extremes[which(Indexes_==1),2]/
  Vect_extremes[which(Indexes_==1),1]
Theta_2<-Vect_extremes[-which(Indexes_==1),1]/
  Vect_extremes[-which(Indexes_==1),2]
## ML gamma
Graphics_estimators_gamma(series = Theta_1,vect_k = c(10:200),
                          Title_graphic = "Extr of Theta_1",
                          NB_years = NULL)
Graphics_estimators_gamma(series = Theta_2,vect_k = c(10:200),
                          Title_graphic = "Extr of Theta_2",
                          NB_years = NULL)

## QQplot

## Highlight the type of distrib+ value of the shape param
qqplot(x=log(Theta_1) ,y=rexp(n = 1000,rate = 1))

List_Params_RF<-list("parametric"=TRUE,"name_RF"=Name_riskF,
                     "function"=Risk_f,"name_model"="log",
                     "Scale_Frechet"=scale_frechet_d,
                     "weightFunction"=wFun,
                     "dWeightFunction"=dwFun,
                     "init_opt_param"=0.5)
## Approach with one or several PCA basis. 
## If one basis is used, choose the number of scores used.
## Otherwise, choose for each forcing cond. the number 
## of scores.
ONE_PCA_BASIS<-TRUE
if(ONE_PCA_BASIS){
  LIST_NB_SCORES<-NULL
  NB_SCORES_OMEGA<-7
}else{
  LIST_NB_SCORES<-list("U"=4,"Surcote"=3)
  NB_SCORES_OMEGA<-NULL
}

# 3) Running simulations --------------------------------------------------
source("functions_MV/functions_Angle_MV.R")
FTF<-function(x){return(log(x))}
FTF_inv<-function(x){return(exp(x))}
re<-Simul_MV_residuals(result_transformation = Result_MV,
                   l_variables = l_name, Q_thresh = Qfound,
                   list_nb_scores = LIST_NB_SCORES,
                   M = M,rotations_available=Bool_rotations,
                   Params_risk_Function=List_Params_RF
                   ,root_for_export=ROOT_export_graph
                   ,One_PCA_base = ONE_PCA_BASIS,
                   NbScores_Omega = NB_SCORES_OMEGA,
                   f_transf=FTF,f_transf_inv=FTF_inv)
U_sim_frech<-re$Frechet_normal$simul[[l_name[[1]]]]
U_obs_frech<-re$Frechet_normal$obs[[l_name[[1]]]]

Obj_sim<-matrix(NA,nrow = nrow(U_sim_frech),
            ncol = ncol(U_sim_frech)*length(l_name))
Obj_data<-matrix(NA,nrow = nrow(U_obs_frech),
                 ncol = ncol(U_obs_frech)*length(l_name))

Root_sim<-re$Frechet_normal$simul
Root_data<-re$Frechet_normal$obs
Z<-0
for(name in l_name){
  Deb<-Z*ncol(U_sim_frech)+1
  END<-(Z+1)*ncol(U_sim_frech)
  Df_sim<-Root_sim[[Z+1]]
  L2_sim<-apply(Df_sim,MARGIN = 1,FUN = calcul_norme_L2)
  Theta_sim<-t(t(Df_sim)%*%diag(L2_sim^(-1)))
  Obj_sim[,c(Deb:END)]<-Theta_sim
  rood<-Root_data[[Z+1]]
  colnames(rood)<-c(1:length(names(rood)))
  Df_data<-as.matrix(rood)
  L2_data<-apply(Df_data,MARGIN = 1,FUN = calcul_norme_L2)
  Theta_data<-t(t(Df_data)%*%diag(L2_data^(-1)))
  Obj_data[,c(Deb:END)]<-Theta_data
  Z<-Z+1
  M_shown<-100
  Filename_omega<-paste0(ROOT_export_graph,"Theta/",name,"/",
                    paste0("Omega_sim_vs_obs_",name,".png"))
  Min_y<-min(c(apply(Theta_data,MARGIN = 2,FUN = min),
               apply(Theta_sim,MARGIN = 2,FUN = min)))
  Max_y<-max(c(apply(Theta_data,MARGIN = 2,FUN = max),
               apply(Theta_sim,MARGIN = 2,FUN = max)))
  png(filename =  Filename_omega,width=1400,
      height=500)
  par(mfrow=c(1,2))
  matplot(t(Theta_data[1:M_shown,]),type="l",
          xlab="Index",ylim=c(Min_y,Max_y),
          ylab="value (-)",cex.lab=1.5)
  matplot(t(Theta_sim[1:M_shown,]),type="l",
          ylim=c(Min_y,Max_y),xlab="Index",
          ylab=NA,cex.lab=1.5)
  par(mfrow=c(1,1))
  dev.off()
}

### Export Theta for future analysis 
write.csv(file=paste0(ROOT_export,"Theta_data.csv"),
          x = Obj_data)
write.csv(file=paste0(ROOT_export,"Theta_sim.csv"),
          x = Obj_sim)
### (09/10) Limitations of PCA for Wind
### U
T<-1
S<-37
# png(filename = paste0(ROOT_export_graph,
#     "Omega_data_sim_t=",T,"_s=",S,".png"),
#     width=1400,height=800)
#par(mfrow=c(1,2))

U<-Obj_data[,T]
V<-Obj_data[,S]
DF_1<-cbind.data.frame(U,V)
DF_sim<-cbind.data.frame(Obj_sim[,T],Obj_sim[,S])
colnames(DF_1)<-sapply(c(1:2),function(x){return(paste0("V",x))})
colnames(DF_sim)<-colnames(DF_1)
df_combo<-rbind.data.frame(DF_1,DF_sim)
df_combo$Legend<-c(rep("data",nrow(DF_1)),
                   rep("simulations",nrow(DF_sim)))
cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")
GG1<-ggplot(data=df_combo,aes(x=V1,y=V2,colour=Legend,shape=Legend),
)+
  geom_point(aes(size=Legend))+
  geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  
  scale_color_manual(values=cols_)+
  scale_fill_manual(values=cols_)+
  scale_shape_manual(values = c("simulations"=17,"data"=19))+
  scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
  guides(fill="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  xlab(paste0(l_name[1],"_",T))+ylab(paste0(l_name[1],"_",S))
ggsave(filename = paste0(Theta_path,
                         paste0("/Distrib_Omega_",l_name[1],"_t=",
                                T,"_s=",S,".png")),plot = GG1,
       width=8,height=6)
#dev.off()

### Surge
talt<-T+37
salt<-S+37
U_data<-Obj_data[,talt]
V<-Obj_data[,salt]
DF_2<-cbind.data.frame(U,V)
DF_sim_2<-cbind.data.frame(Obj_sim[,talt],Obj_sim[,salt])
colnames(DF_2)<-sapply(c(1:2),function(x){return(paste0("V",x))})
colnames(DF_sim_2)<-colnames(DF_2)
df_combo2<-rbind.data.frame(DF_2,DF_sim_2)
df_combo2$Legend<-c(rep("data",nrow(DF_2)),
                   rep("simulations",nrow(DF_sim_2)))
GG2<-ggplot(data=df_combo2,aes(x=V1,y=V2,colour=Legend,shape=Legend),
)+
  geom_point(aes(size=Legend))+
  geom_xsidedensity(data=df_combo2,aes(fill=Legend), alpha = 0.5)+
  geom_ysidedensity(data=df_combo2,aes(fill=Legend), alpha = 0.5)+
  scale_color_manual(values=cols_)+
  scale_fill_manual(values=cols_)+
  scale_shape_manual(values = c("simulations"=17,"data"=19))+
  scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
  guides(fill="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  xlab(paste0(l_name[2],"_",talt))+ylab(paste0(l_name[2],"_",salt))
ggsave(filename = paste0(Theta_path,
                         paste0("/Distrib_Omega_",l_name[2],"_t=",
                                talt,"_s=",salt,".png")),plot = GG2,
       width=8,height=6)

# plot(U,V)
# points(Obj_sim[,talt],Obj_sim[,salt],col="red"
#        ,cex=0.5)
# par(mfrow=c(1,1))
# dev.off()
### Atom on 0 for U ?? 

mu_Theta<-apply(X = Obj_data,MARGIN = 2,FUN = mean)
N_rep<-500
B<-nrow(Obj_data)
N_rep<-500
lq<-c(0.05,0.50,0.95,0.975)
Tend_Theta<-replicate(n =N_rep,expr = Resamples_tendencies(B = B,
                                                     extremes_indus=Obj_data,
                                                     list_Q =lq))
Results_data<-sapply(X = lq,FUN = function(x,q){return(apply(x,MARGIN = 2,FUN =
                                                 function(x){return(quantile(x,q))}))},
       x=Obj_data)
Results_sim<-sapply(X = lq,FUN = function(x,q){return(apply(x,MARGIN = 2,FUN =
                                                               function(x){return(quantile(x,q))}))},
                     x=Obj_sim)
bound_minus05<-apply(Tend_Theta[1,,],MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.025))})
bound_plus05<-apply(Tend_Theta[1,,],MARGIN = 1,
                    FUN = function(x){return(quantile(x,0.975))})

##########
bound_minus50<-apply(Tend_Theta[2,,],MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.025))})
bound_plus50<-apply(Tend_Theta[2,,],MARGIN = 1,
                    FUN = function(x){return(quantile(x,0.975))})

########
bound_minus95<-apply(Tend_Theta[3,,],MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.025))})
bound_plus95<-apply(Tend_Theta[3,,],MARGIN = 1,
                    FUN = function(x){return(quantile(x,0.975))})

######
bound_minus975<-apply(Tend_Theta[4,,],MARGIN = 1,
                      FUN = function(x){return(quantile(x,0.025))})
bound_plus975<-apply(Tend_Theta[4,,],MARGIN = 1,
                     FUN = function(x){return(quantile(x,0.975))})

# Graphics ----------------------------------------------------------------

df<-cbind.data.frame(c(Results_data),
                     c(bound_minus05,bound_minus50,
                       bound_minus95,bound_minus975),
                     c(bound_plus05,bound_plus50,
                       bound_plus95,bound_plus975))
colnames(df)<-c("estimator","bound_minus","bound_plus")

df$estimator_simul<-c(Results_sim)
df$Time<-rep(c(1:ncol(Obj_data)),4)
df$percent<-c(rep("Q05",ncol(Obj_data)),
              rep("Q50",ncol(Obj_data)),
              rep("Q95",ncol(Obj_data)), 
              rep("Q975",ncol(Obj_data)))
df$Time
First<-which(df$Time<=37)
df_first<-df[First,]
GG_percent<-ggplot(data=df_first,aes(x=Time,y=estimator,group=interaction(percent),col="data"))+
  facet_wrap(~percent,scales="free_y")+
  geom_line()+
  geom_point()+
  geom_ribbon(mapping = aes(ymin=bound_minus,ymax=bound_plus,col="confidence_band"),alpha=0.15,
              fill="grey", linetype = "dashed")+
  geom_line(aes(x=Time,y=estimator_simul,col="simulations"),linetype=2)+
  theme(axis.title=element_text(size=25),
        legend.text=element_text(size=14),
        legend.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12))+
  scale_color_manual(values=cols_)+
  xlab("Time")+
  labs(col="Legend")
ggsave(plot = GG_percent,
       filename = paste0(Theta_path,"/",l_name[1],"/Percent_theta_",l_name[1],
                         ".png"),
       width=8,height = 6)

Second<-which(df$Time>37)
df_second<-df[Second,]
GG_percent<-ggplot(data=df_second,aes(x=Time,y=estimator,group=interaction(percent),col="data"))+
  facet_wrap(~percent,scales="free_y")+
  geom_line()+
  geom_point()+
  geom_ribbon(mapping = aes(ymin=bound_minus,ymax=bound_plus,col="confidence_band"),alpha=0.15,
              fill="grey", linetype = "dashed")+
  geom_line(aes(x=Time,y=estimator_simul,col="simulations"),linetype=2)+
  theme(axis.title=element_text(size=25),
        legend.text=element_text(size=14),
        legend.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12))+
  scale_color_manual(values=cols_)+
  xlab("Time")+
  labs(col="Legend")
ggsave(plot = GG_percent,
       filename =paste0(Theta_path,"/",l_name[2],"/Percent_theta_",l_name[2],
                        ".png"),width=8,height = 6)
dev.off()
#dev.off()
j_coord<-1
j_nl2<-1
pval_indep<-cor.test(re$coords_data[,j_coord],
         y = Vect_l_function[re$Indices_exts,j_nl2],
         method = "kendall")$p.value
pval_indep
Sub<-cbind(re$coords_data[,1],re$coords_data[,2])
UNifs<-VineCopula::pobs(Sub)

plot(UNifs[,1],UNifs[,2])

Result<-lm(UNifs[,2]~UNifs[,1]+I(UNifs[,1]^2))
png("test.png")
plot(UNifs[,1],UNifs[,2])
points(UNifs[,1],Result$fitted.values,col="red")
dev.off()

# Theta_found vs cost evolution ---------------------------------
######
if(List_Params_RF$name_model=="log"){
  Chi_estim<-2-2^(re$Param_found)
  print(paste0("The Chi would be around ", round(Chi_estim,2)))
  png(paste0(Lg_path,"/log_model_param_vs_data.png"),
      height=800,width = 1000)
  par(mfrow=c(1,2))
  POT::chimeas(data = Result_MV$df,ask = FALSE,which = 1)
  abline(h=Chi_estim,col="red")
  abline(h=0,col="black")
  POT::chimeas(data = Result_MV$df,ask = FALSE,which = 2)
  par(mfrow=c(1,1))
  dev.off()
}
if(List_Params_RF[["name_RF"]]!="max"){
  Excedents_lprime<-Vect_l_function[re$Indices_exts,]
  Evol_cost<-function(theta_candidate){
    Objective_function(theta = theta_candidate,
                       weightFun = List_Params_RF[["weightFunction"]],
                       dWeightFun = List_Params_RF[["dWeightFunction"]],
                       Extreme_inds =as.list(as.data.frame(t(Excedents_lprime))) ,
                       u = re$threshold,
                       model_used = List_Params_RF[["name_model"]])

  }
  vect_candidates<-seq.int(from = 0.6,to = 0.99,length.out = 20 )
  Cost_evol<-sapply(vect_candidates,
                    Evol_cost)
  # Export evolution alpha --------------------------------------------------
  png(file = paste0(ROOT_export,"evol_alpha_param.png"),
      bg = "transparent")
  plot(vect_candidates,Cost_evol,xlab="Param",ylab="Cost")
  dev.off()
}

# Summary of the results --------------------------------------------------
Riskmv<-re$Risk_MV
Inds_exts<-re$Indices_exts
### if you use EXTGPD, you must only take pos values-->
Inds_exts<-Inds_dates_pos[Inds_exts]
DF_inds_risk<-cbind.data.frame(Riskmv,
                               Inds_exts)
colnames(DF_inds_risk)<-c("RiskF","indexes")
#Export dates of extreme events ------------------------------------------
write.csv(file=paste0(ROOT_export,"ext_series_chosen_risks.csv"),
          x = DF_inds_risk)

L_extremes<-Vect_l_function[re$Indices_exts,]
# 
# Kokozka -----------------------------------------------------------------
# #######
source("functions_MV/Fonctions_consistency.R")
Extreme_corr_data<-sapply(c(1:length(re$Indices_exts)),Extreme_corr,
       liste_MV_simul=re$Frechet_normal$obs,l_name=l_name,
       L=37,d=2)
N<-nrow(re$Frechet_normal$simul[[ind_var1]])
Extreme_corr_simul<-sapply(c(1:N),Extreme_corr,
                          liste_MV_simul=re$Frechet_normal$simul,l_name=l_name,
                          L=37,d=2)
Mu_RHO<-mean(Extreme_corr_data)
sd_RHO<-sd(Extreme_corr_data)
print(c(Mu_RHO,sd_RHO))
Mu_simul_rho<-mean(Extreme_corr_simul)
sd_simul_rho<-sd(Extreme_corr_simul)
print(c(Mu_simul_rho,sd_simul_rho))

Maxi<-apply(X = Result_MV$df,MARGIN = 1,FUN = max)
inds_exts<-which(Maxi>re$threshold)
z<-1
marg_z<-Result_MV$df[inds_exts,z]

# # # l function --------------------------------------------------------------
# Law of the couple ------------------------------------------
# j<-1
# Marg<-as.numeric(sapply(Vect_l_function[re$Indices_exts,j],CDF_1d_Mbiv,Th = re$threshold,
#                         Params_Biv = Mod_BIV,
#                         z = j))
# goftest::ad.test(x = Marg,null = "punif")$p.value

# Shapes and coordinates --------------------------------------------------
par(mfrow=c(1,2))
#ind_var1

if(ONE_PCA_BASIS==FALSE){
  MAX_Eigen<-max(re$EIGEN_functions[[ind_var1]])
  MIN_Eigen<-min(re$EIGEN_functions[[ind_var1]])
  plot(c(1:37),re$EIGEN_functions[[ind_var1]][,1],ylim=c(MIN_Eigen,MAX_Eigen),
       ylab=expression(nu[j]),xlab="Time")
  lines(c(1:37),re$EIGEN_functions[[ind_var1]][,2],col="green")
  lines(c(1:37),re$EIGEN_functions[[ind_var1]][,3],col="red")
  
  MAX_Eigen<-max(re$EIGEN_functions[[ind_var2]])
  MIN_Eigen<-min(re$EIGEN_functions[[ind_var2]])
  plot(c(1:37),re$EIGEN_functions[[ind_var2]][,1],ylim=c(MIN_Eigen,MAX_Eigen),
       ylab=expression(nu[j]),xlab="Time")
  lines(c(1:37),re$EIGEN_functions[[ind_var2]][,2],col="green")
  mtext("PCA eigen functions for the variables",outer=TRUE,
        line=-2)
  
}else{
  MAX_Eigen<-max(re$EIGEN_functions)
  MIN_Eigen<-min(re$EIGEN_functions)
  plot(c(1:37),re$EIGEN_functions[c(1:37),1],ylim=c(MIN_Eigen,MAX_Eigen),
       ylab=expression(nu[j]),xlab="Time")
  lines(c(1:37),re$EIGEN_functions[c(1:37),2],col="green")
  lines(c(1:37),re$EIGEN_functions[c(1:37),3],col="red")

  plot(c(1:37),re$EIGEN_functions[c(38:74),1],ylim=c(MIN_Eigen,MAX_Eigen),
       ylab=expression(nu[j]),xlab="Time")
  lines(c(1:37),re$EIGEN_functions[c(38:74),2],col="green")
  lines(c(1:37),re$EIGEN_functions[c(1:37),3],col="red")
}
png(filename = paste0(Theta_path,"/corrs_BW_coords_sim_data.png"),
    height=800,width=1000)
par(mfrow=c(1,2))
corrplot::corrplot(cor(re$coords_data))
corrplot::corrplot(cor(re$coords_simul))
par(mfrow=c(1,1))
dev.off()
j<-5
plot(density(re$coords_data[,j]))
plot(density(re$coords_simul[,j]))
par(mfrow=c(1,1))

# Coords ------------------------------------------------------------------
NDIM<-ncol(re$coords_simul)
KS_results_coords<-rep(NA,NDIM)
for(j in c(1:NDIM)){
  KS_results_coords[j]<-ks.test(re$coords_simul[,j],
                                re$coords_data[,j])$p.value
}
write.csv(file = paste0(ROOT_export,"coord_data.csv"),
          x = re$coords_data)
write.csv(file = paste0(ROOT_export,"coord_simul.csv"),
          x = re$coords_simul)
KS_results_coords
# Illustration ------------------------------------------------------------
Simul<-re$simul

N_shown<-10
Int_conf<-function(df,Q,df_simul){
  int_1<-apply(X = df,MARGIN=2,FUN = function(x){
    return(quantile(x,Q))
  })
  int_2<-apply(X = df_simul,MARGIN=2,FUN = function(x){
    return(quantile(x,Q))
  })
  return(list("obs"=int_1,
              "simul"=int_2))
}
Int_surge<-Int_conf(df = Simul[[ind_var1]],
                    Q = 0.025,df_simul = re$obs_exts[[ind_var1]])
Int_surge2<-Int_conf(df = Simul[[ind_var1]],
                    Q = 0.975,df_simul = re$obs_exts[[ind_var1]])
Int_Hs<-Int_conf(df = Simul[[ind_var2]],
                    Q = 0.025,df_simul = re$obs_exts[[ind_var2]])
Int_Hs2<-Int_conf(df = Simul[[ind_var2]],
                     Q = 0.975,df_simul = re$obs_exts[[ind_var2]])
DF_first<-rbind.data.frame(Simul[[ind_var1]],
                           re$obs_exts[[ind_var1]])

Indexes_data<-sample(c(1:nrow(Simul[[ind_var1]])),size = N_shown,
                     replace=FALSE)
Indexes_sim<-sample(c(1:nrow(re$obs_exts[[ind_var1]])),size = N_shown,
                     replace=FALSE)

Source<-c(rep("simul",nrow(Simul[[ind_var1]])),
              rep("data",length(re$Indices_exts)))
DF_first$source<-Source
DF_second<-rbind.data.frame(Simul[[ind_var2]],re$obs_exts[[ind_var2]])
DF_second$source<-Source

Melteur_first<-melt(DF_first)
L<-table(Melteur_first$variable)[1]
Melteur_first$index<-as.character(rep(c(1:L),nrow(Melteur_first)/L))
Melteur_first$variable<-as.numeric(Melteur_first$variable)
colnames(Melteur_first)<-c("origin","Time",
                            "value","index")
Melteur_first$ymin<-ifelse(Melteur_first$origin=="simul",
                           Int_surge$simul[Melteur_first$Time],
                           Int_surge$obs[Melteur_first$Time])
Melteur_first$ymax<-ifelse(Melteur_first$origin=="simul",
                           Int_surge2$simul[Melteur_first$Time],
                           Int_surge2$obs[Melteur_first$Time])
Melteur_second<-melt(DF_second)

L<-table(Melteur_second$variable)[1]
Melteur_second$index<-as.character(rep(c(1:L),nrow(Melteur_second)/L))
Melteur_second$variable<-as.numeric(Melteur_second$variable)
colnames(Melteur_second)<-c("origin","Time",
                            "value","index")
Melteur_second$ymin<-ifelse(Melteur_second$origin=="simul",
                           Int_Hs$simul[Melteur_first$Time],
                           Int_Hs$obs[Melteur_first$Time])
Melteur_second$ymax<-ifelse(Melteur_second$origin=="simul",
                           Int_Hs2$simul[Melteur_first$Time],
                           Int_Hs2$obs[Melteur_first$Time])

Melteur_all<-rbind.data.frame(Melteur_first,Melteur_second)
Melteur_all$functional<-c(rep(l_name[1],nrow(Melteur_first)),
                        rep(l_name[2],nrow(Melteur_second)))


D_<-which(Melteur_all$origin=="data")
Ind_used_data<-sample(Melteur_all[D_,]$index,size = N_shown,
                      replace = FALSE)
Ind_used_sim<-sample(Melteur_all[-D_,]$index,size = N_shown,
                     replace = FALSE)

cond_data<-Melteur_all$index%in%union(Ind_used_sim,
                                      Ind_used_data)
Melteur_all$Time<-(Melteur_all$Time-19)/6
Melteur_all$index<-as.numeric(Melteur_all$index)
EXPLES_simul<-ggplot(Melteur_all[cond_data,],aes(x=Time,y=value,col=index,
                         group=interaction(index)))+
  guides(col="none")+
  facet_wrap(~functional+origin)+
  geom_line()+
  scale_color_continuous(type="viridis")
EXPLES_simul
ggsave(filename=paste0(TEND_path,"/exples_simul_MV.png"),
       plot = EXPLES_simul,width = 6,height = 6)

Mlt_copy<-Melteur_all
Mlt_copy$Time<-as.character(Mlt_copy$Time)

# Violin plot -------------------------------------------------------------
############
Mlt_copy$Time<-as.character(Mlt_copy$Time)
DISTRIB_marginals_plot<-ggplot(Mlt_copy,aes(y=value,x=Time,fill=functional))+
  geom_boxplot(alpha=0.2,width=0.1)+
  geom_violin()+
  xlab(" ")+
  facet_wrap(~functional+origin,scales="free_y")+
  labs(fill="Legend")+
  theme(axis.text.x = element_blank())
  
ggsave(filename=paste0(ROOT_export_graph,"Tend/violin_plot.png"),
       plot = DISTRIB_marginals_plot,width = 6,
       height = 6)

par(mfrow=c(1,1))
Pareto_S<-Result_MV$resume[[ind_var1]]$transf
Pareto_Hs<-Result_MV$resume[[ind_var2]]$transf
L2_P_S<-apply(Pareto_S[re$Indices_exts,],MARGIN = 1,FUN = calcul_norme_L2)
L2_P_HS<-apply(Pareto_Hs[re$Indices_exts,],MARGIN = 1,FUN = calcul_norme_L2)
ks.test(L2_P_HS,L2_P_S)
wilcox.test(L2_P_HS,L2_P_S)

L2_pto_sur<-apply(re$Frechet_normal$simul[[ind_var1]],MARGIN = 1,
                  FUN = calcul_norme_L2)
L2_pto_Hs<-apply(re$Frechet_normal$simul[[ind_var2]],MARGIN = 1,
                 FUN = calcul_norme_L2)
DF_data_vectL<-matrix(NA,ncol=length(names(re$Frechet_normal$simul)),
                      nrow =length(re$Indices_exts))
for(j in c(1:length(l_name))){
  name_j<-names(Result_MV$resume)[j]
  PTO_j<-Result_MV$resume[[name_j]]$transf
  Output_j<-apply( PTO_j[re$Indices_exts,],MARGIN = 1,
                            FUN = calcul_norme_L2)
  DF_data_vectL[,j]<-Output_j
}
# Z<-2
# P_BIV<-re$modele_maj
# R_X<-sapply(L2_pto_Hs,CDF_1d_Mbiv,z=Z,
#             Params_Biv=P_BIV,Th=re$threshold)
# summary(R_X)
# Index<-which(is.na(R_X)==TRUE)
# Subset<-R_X
# if(length(Index)>0){
#   Subset<-R_X[-Index]
# }
# plot(density(Subset))
# goftest::ad.test(Subset,null = "punif")

Min_x<-min(min(L2_pto_sur),min(L2_P_S))
Max_x<-max(max(L2_pto_sur),max(L2_P_S))

Min_y<-min(min(L2_pto_Hs),min(L2_P_HS))
Max_y<-max(max(L2_pto_Hs),max(L2_P_HS))
X_name<-ifelse(l_name[1]=="Surcote",yes = "Surge",no=l_name[1])
Y_name<-ifelse(l_name[2]=="Surcote",yes = "Surge",no=l_name[2])
par(mfrow=c(1,1))
png(filename = paste0(Lg_path,"/pair_Lg.png"))

plot(L2_P_S,L2_P_HS,log = 'xy',
     xlab=paste0("Z(",X_name,") "),
     ylab=paste0("Z(",Y_name,") "),
     ylim = c(Min_y,Max_y),
     xlim=c(Min_x,Max_x))
points(L2_pto_sur,
       L2_pto_Hs,col="red",cex=0.5)
if(!Name_riskF%in%c("Weighted_NQ","max","min")){
  ## Change the threshold used.
  J_plot<-2
  X_private<-DF_data_vectL[,-J_plot]
  x_bound<-seq.int(from = 0.001,to =20,length.out =100)
  y_bound<-sapply(x_bound,FUN=Shape_boundary,
                 j = J_plot,ul = re$threshold,RiskF = Risk_f)
  lines(x_bound,y_bound,col="blue")
}else{
  abline(h=re$threshold,v=re$threshold,col="blue")
}
dev.off()

## originally
Min_x2<-min(min(Vect_l_function[,1]),min(L2_pto_sur))
Max_x2<-max(max(Vect_l_function[,1]),max(L2_pto_sur))

Min_y2<-min(min(Vect_l_function[,2]),min(L2_pto_Hs))
Max_y2<-max(max(Vect_l_function[,2]),max(L2_pto_Hs))
plot(Vect_l_function,xlim=c(Min_x2,Max_x2),
     ylim=c(Min_y2,Max_y2),log="xy")

DATES<-read.csv(file = "../ss_tend/HIVER/dates_prises.csv")[,2]

Inds_2008<-which(DATES=="10/03/2008")
Dates_storm<-c("10/03/2008","28/02/2010","10/02/2009",
               "27/10/2004","07/02/2001","10/01/2001")
Cplmt_dates_1<-c("26/02/2010","08/03/2008","08/02/2009",
              "25/10/2004","05/02/2001","08/01/2001")

Cplmt_dates_2<-c("09/03/2008","27/02/2010","09/02/2009",
                 "26/10/2004","06/02/2001","09/01/2001")
Inds_storm<-union(which(substr(DATES,1,10)%in%Dates_storm),
                  which(substr(DATES,1,10)%in%Cplmt_dates_1)
                  ) 
Inds_storm<-union(Inds_storm,
                  which(substr(DATES,1,10)%in%Cplmt_dates_2)
                  )
DATES[Inds_storm]
ind_Joh<-which(substr(DATES,1,10)=="10/03/2008")
DATA_pos_neg<-read.csv(file = paste0(l_root,l_name[1],"_residuals.csv"))[,c(2:38)]
DATA_pos_neg[ind_Joh,]
### select storm found in the extreme events seen
ind_storm_found<-intersect(Inds_storm,
                           Result_MV$SUB_set)
Storms_seen<-DATES[ind_storm_found]
F_return_index<-function(x){
  return(Result_MV$SUB_set[x])
}
Index_storm<-which( Result_MV$SUB_set%in%ind_storm_found)
Correspondg_vectors<-Vect_l_function[Index_storm,]
df_sim<-cbind.data.frame(L2_pto_sur,L2_pto_Hs)
colnames(df_sim)<-sapply(c(1:length(l_name)),
                         FUN=function(x){
                           return(paste0("V",x))
                         })
VECT2<-Vect_l_function
colnames(VECT2)<-colnames(df_sim)
Correspondg_vectors
X_storms<-Correspondg_vectors[,1]
Y_storms<-Correspondg_vectors[,2]
ANOTATE_data<-data.frame(x=X_storms,y=Y_storms,
                         label=rep("storm",length(Index_storm)))

GG_data_and_exceed<-ggplot(data=VECT2,aes(x=V1,y=V2))+
  geom_point()+
  geom_point(data = ANOTATE_data, 
             aes(x, y,col="storm"), size = 4,
             shape=0) +
  annotate("text", x = (X_storms)*(8/10), y = Y_storms*(9/10), 
           label=Storms_seen,colour = "red", 
           size = 5)+
  scale_y_continuous(transform="log",
                     labels = label_number(accuracy = 0.1))+
  scale_x_continuous(transform="log",
                     labels = label_number(accuracy = 0.1))+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=14),
        legend.title = element_text(size=15),
        axis.text = element_text(size=14))+
  geom_hline(aes(yintercept = re$threshold,col="exceedances"),
             linetype="dashed")+
  geom_vline(aes(xintercept=re$threshold,col="exceedances"),
             linetype="dashed")+
  scale_color_manual(values=c("exceedances"="red",
                              "storm"="black"))+
  labs(col="Legend")+
  xlab(paste0("Norm of T(",l_name[1],") (.)"))+
  ylab(paste0("Norm of T(",l_name[2],") (.)"))

ggsave(plot = GG_data_and_exceed,
       filename = paste0(Lg_path,"/norm_ext.png"),
       width = 8,height = 8)
DF_sim_vs_obs<-rbind.data.frame(VECT2,df_sim)
DF_sim_vs_obs$origin<-c(rep("data",nrow(VECT2)),
                        rep("simulations",nrow(df_sim)))

GG_repre_exts_norms<-ggplot(data=DF_sim_vs_obs,aes(x=V1,y=V2))+
  geom_point(aes(col=origin,size=origin))+
  geom_point(data = ANOTATE_data, 
             aes(x, y,col=label),size=4,shape=0) +
  annotate("text", x = (X_storms)*(8/10), y = Y_storms*(9/10), 
           label=Storms_seen,colour = "red", 
           size = 5)+
  scale_y_continuous(transform="log",
                     labels = label_number(accuracy = 0.1))+
  scale_x_continuous(transform="log",
                     labels = label_number(accuracy = 0.1))+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=14),
        legend.title = element_text(size=15),
        axis.text = element_text(size=14))+
  geom_hline(aes(yintercept = re$threshold,col="exceedances"),
             linetype="dashed",)+
  geom_vline(aes(xintercept=re$threshold,col="exceedances"),
             linetype="dashed")+
  scale_color_manual(values=c("data"="blue",
                              "simulations"="orange",
                              "storm"="black",
                              "exceedances"="red"))+
  scale_size_manual(values=c("data"=2,
                             "simulations"=1))+
  # scale_shape_manual(values=c("data"=2,
  #                            "simulations"=2,
  #                            "storm"=0))+
  labs(col="Legend")+
  guides(size="none")+
  xlab(paste0("Norm of T(",l_name[1],") (.)"))+
  ylab(paste0("Norm T(",l_name[2],") (.)"))

ggsave(plot = GG_repre_exts_norms,
       filename = paste0(Lg_path,"/ext_simuls_RM.png"),
       width = 8,height = 8)

# legend(x="topright",legend=c("data","simulations"),
#        col=c(1,2),lwd=2,bty="n")
cor(L2_P_S,L2_P_HS,method = "kendall")
cor(L2_pto_sur,L2_pto_Hs
    ,method = "kendall")
# Melt_couple_l<-rbind.data.frame(L_ext_d,L_simul)
# require(scales)
# ggplot(data=Melt_couple_l,aes(x=Surcote,y=
#                                 Hs,colour=Legend,shape=Legend))+
#   geom_point(aes(size=Legend))+
#   geom_xsidedensity(data=Melt_couple_l,aes(fill=Legend), alpha = 0.5)+
#   geom_ysidedensity(data=Melt_couple_l,aes(fill=Legend), alpha = 0.5)+
#   scale_color_manual(values=cols_)+
#   scale_fill_manual(values=cols_)+
#   scale_shape_manual(values = c("simulations"=17,"data"=19))+
#   scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
#   scale_y_continuous(labels = custom_labels,transform = "log")+
#   scale_x_continuous(labels = custom_labels,transform = "log")+
#   guides(fill="none")

POT::chimeas(cbind(L2_pto_Hs,L2_pto_sur),ask=FALSE)

Maximum<-apply(cbind(L2_pto_Hs,L2_pto_sur),MARGIN = 1,FUN = max)
summary(Maximum)

Maximum_obs<-apply(cbind(L2_P_S,L2_P_HS),MARGIN = 1,FUN = max)
summary(Maximum_obs)

# MV KS TEST -------------------------------------------------------------------
# Kernel_dens_d<-function(vect_2d,h,new_x_y){
#   Values<-(apply((vect_2d-new_x_y)**2,MARGIN = 1,FUN = sum))**(1/2)
#   Values_dens<-sapply(Values,FUN = function(x){
#     return(dnorm(x,sd = h))
#   })
#   return(mean(Values_dens))
# }
df<-cbind(L2_P_HS,L2_P_S)

# f_hat_values<-apply(df,FUN=Kernel_dens_d,vect_2d =df,
#               h = length(L2_P_HS)^(-4/3),
#               MARGIN = 1)
M<-2000
Unif_simul<-runif(M)
z<-1
# Sim_l<-t(parSapply(cl = coeurs,Unif_simul,Function_one_couple_l,
#                 Params_Biv = P_BIV,z = z,
#                 Th = re$threshold))
# #with data -> df.
# F_theo_values<-apply(Sim_l,FUN =CDF_law_2d,
#                      Th=re$threshold,
#                      Params_Biv =P_BIV,
#                      MARGIN = 1)
# summary(F_theo_values)
# summary(Estimation_values)

fUl<-min(apply(cbind(L2_P_HS,L2_P_S),MARGIN = 1,FUN = max))
cor(L2_P_HS,L2_P_S,method = "kendall")
cor(L2_pto_Hs,L2_pto_sur,method = "kendall")

# PCA for the angle -------------------------------------------------------
val_chosen<-c(3,2)
NPCA<-as.list(val_chosen)


#Calculate correlations --------------------------------------------------
#Shape correlation var1 var2 ---------------------------------------------
Mat_obs<-Fct_correlations(method_corr = "kendall",
                        df1 =Obj_data[,c(1:37)],
                        df2=Obj_data[,c(38:74)])
apply(Obj_data[,c(1:37)],MARGIN = 2,FUN = sd)
apply(Obj_sim[,c(1:37)],MARGIN = 2,FUN = sd)
Conf_interval_angle<-Conf_interval_correlations_2V(df1 = Obj_data[,c(1:37)],
                                                  df2 = Obj_data[,c(38:74)],
                                                  niv_conf = 0.05,
                                                  method_corr = "kendall",
                                                  NB_boot = 500)
Mat_simul<-Fct_correlations(method_corr = "kendall",
                            df1 =Obj_sim[,c(1:37)],
                            df2=Obj_sim[,c(38:74)])
Data_corr<-as.data.frame(Conf_interval_angle)
colnames(Data_corr)<-c("bound_inf","data","bound_sup")
Data_corr$simulations<-Mat_simul
Data_corr$Time<-(c(1:nrow(Data_corr))-19)/6
GG_corr_shape<-ggplot(data=Data_corr,aes(x=Time,y=data,col="data"))+
  geom_ribbon(aes(ymin=bound_inf,col="confidence_band",
                  ymax=bound_sup),alpha=0.05,fill="grey",
              linetype="dashed")+
  geom_line(aes(y=simulations,col="simulations"))+
  geom_line()+
  geom_point()+
  xlab("Time (hour) with respect to tidal peak")+
  scale_color_manual(values = cols_)+
  ylab("correlation")+
  labs(col="Legend")+theme(axis.title=element_text(size=15),
                            legend.text=element_text(size=10))
ggsave(filename = paste0(ROOT_export_graph,"corr_shape.png"),
       plot = GG_corr_shape,width=6,height = 6)
#dev.off()
# Correlations ------------------------------------------------------------
Surcote_Frechet<-Result_MV$resume[[ind_var1]]$transf[re$Indices_exts,]
Hs_Frechet<-Result_MV$resume[[ind_var2]]$transf[re$Indices_exts,]
Corr_Frechet<-Fct_correlations(method_corr = "kendall",
                               df1 =Surcote_Frechet,
                               df2=Hs_Frechet)
Conf_interval_data_Frechet<-Conf_interval_correlations_2V(
                              df1 = Surcote_Frechet,
                              df2 = Hs_Frechet,
                              niv_conf = 0.05,
                              method_corr = "kendall",
                              NB_boot = 500)
Hs_RV<-re$Frechet_normal$simul[[ind_var2]]
Sur_RV<-re$Frechet_normal$simul[[ind_var1]]
Corr_sim_Frechet<-Fct_correlations(method_corr = "kendall",
                                   df1 =Sur_RV,
                                   df2=Hs_RV)
Data_corr_whole<-as.data.frame(Conf_interval_data_Frechet)
colnames(Data_corr_whole)<-c("bound_inf","data","bound_sup")
Data_corr_whole$simulations<-Corr_sim_Frechet
Data_corr_whole$Time<-(c(1:nrow(Data_corr_whole))-19)/6

Corr_ggplot<-ggplot(data=Data_corr_whole,aes(x=Time,y=data,col="data"))+
  geom_ribbon(aes(ymin=bound_inf,col="confidence_band",
                  ymax=bound_sup),alpha=0.05,fill="grey",
              linetype="dashed")+
  geom_line(aes(y=simulations,col="simulations"))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cols_)+
  xlab("Time (hour) with respect to tidal peak")+
  ylab("correlation")+
  labs(col="Legend")+theme(axis.title=element_text(size=15),
                           legend.text=element_text(size=10))
ggsave(filename = paste0(ROOT_export_graph,"corr_Kendall.png"),
       plot = Corr_ggplot,width=6,height = 6)

# (Shapes) Correlation between variables OK -------------------------------
######################
liste_ext<-Create_list_fromALL(l_name = l_name,
                obj_source = re$obs_exts)
# Correlations Pareto scale -----------------------------------------------
Mat_simul_Pareto<-Fct_correlations(method_corr = "kendall",
                            df1 =re$simul[[ind_var1]],
                            df2=re$simul[[ind_var2]])

Mat_obs_Pareto<-Fct_correlations(method_corr = "kendall",
                          df1 =re$obs_exts[[ind_var1]],
                          df2=re$obs_exts[[ind_var2]])
summary(Mat_simul_Pareto)
summary(Mat_obs_Pareto)
plot(c(1:37),y = Mat_obs_Pareto,ylim=c(-0.5,0))
lines(c(1:37),y = Mat_simul_Pareto,col="red")
legend(x="topright",legend=c("data","simul"),
       col=c(1,2),lwd=2,bty="n")

# Correlation var1 var2 ---------------------------------------------------

Mat_corr_simul<-Fct_correlations(method_corr = "kendall",
                           df1 =re$simul[[ind_var2]],
                           df2=re$simul[[ind_var1]])
Mat_corr<-Fct_correlations(method_corr = "kendall",
                           df1 =liste_ext[[ind_var2]],
                           df2=liste_ext[[ind_var1]])
plot(c(1:37),y = Mat_corr,
     main="Correlations between Hs and S at each time (red in simul)")
lines(c(1:37),y = Mat_corr_simul,
      col="red")


# Marginal law ------------------------------------------------------------
NB_times<-ncol(re$simul[[ind_var1]])
KS_test_S<-rep(NA,NB_times)
KS_test_HS<-KS_test_S
for(t in c(1:NB_times)){
  KS_test_S[t]<-ks.test(x =re$simul[[ind_var1]][,t],
                      y=liste_ext[[ind_var1]][,t])$p.value
  KS_test_HS[t]<-ks.test(x = re$simul[[ind_var2]][,t],
                         y=liste_ext[[ind_var2]][,t])$p.value
}
plot(KS_test_S)
plot(KS_test_HS)
### on conserve l'hypothese nulle ! ###
# cross extremogram -------------------------------------------------------
# Obj<-Extremo_Quality_MV(list_simul = re$simul,
#                    list_reality = liste_ext,
#                    name_other = l_name[2],
#                    name_cond = l_name[1],
#                    q_per_time = 0.95,
#                    B = 500)
# 
# ggplot(data = Obj,aes(x=time_d,y=val_data,col="data"))+
#   geom_line()+
#   geom_point()+
#   geom_line(aes(y=val_sim,col="simulations"))+
#   geom_ribbon(aes(ymin=Qinf,ymax=Qsup,col="confidence_band"),
#               fill="grey",alpha=0.05,linetype="dashed")+
#   ylab("value")+xlab("Lag h")+
#   scale_color_manual(values=cols_)+
#   theme(axis.title=element_text(size=15),
#         legend.text=element_text(size=10))+
#   labs(col="Legend")
  
# Export ------------------------------------------------------------------
##################
for(nameV in l_name){
  write.csv(x =re$simul[[nameV]],file=paste0(ROOT_export,nameV,"_simul.csv"))
  write.csv(x =liste_ext[[nameV]],file=paste0(ROOT_export,nameV,"_obs_ext.csv"))
}

#Fin code ----------------------------------------------------------------
###################
stopCluster(coeurs)

