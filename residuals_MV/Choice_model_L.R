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

### Source des fonctions. ####
source("../fonctions/fonctions_Pareto.R")
source("../fonctions/fonctions.r")
source("../fonctions/fonctions_r_pareto_EVOLTAU.r")
source("functions_MV/function_simul_MV.r")
source("functions_MV/function_estim_simul_gPareto.R")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

# Lancement fonction ------------------------------------------------------
# ------------------------------------------------------------------------

# for(j in c(1:ncol(TS_source))){
#   sigma_egpd_t<-Model_egpd$sigma
#   Gamma_egpd_t<-Model_egpd$mu
#   Nu_t<-Model_egpd$nu
#   db<-as.data.frame(TS_source[,j])
#   colnames(db)<-c("x")
#   con <- gamlss.control(n.cyc = 100,
#                       mu.step = 0.1, 
#                       sigma.step = 0.1, nu.step = 0.1,tau.step = 0.1,autostep=TRUE)
#   con.i=glim.control(glm.trace = FALSE)
#   Fitting_GAMMA_t <- gamlss(x~1, 
#                          data=db,
#                          family = GG(),
#                          control = con,mu.start=exp(0),
#                          sigma.start=exp(sigma_egpd_t),
#                          nu.start=as.numeric(Gamma_egpd_t),
#                          i.control=con.i,
#                          method=CG())
#   plot(Fitting_GAMMA_t)
# }
# print(j)
# MUGAMMA<-as.numeric(Fitting_GAMMA_t$mu.fv)[1]
# NUGAMMA<-as.numeric(Fitting_GAMMA_t$nu.fv)[1]
# SIGMAGAMMA<-as.numeric(Fitting_GAMMA_t$sigma.fv)[1]
# Result_found<-pGG(q = TS_source[,j],mu = MUGAMMA,
#             sigma = SIGMAGAMMA,nu = NUGAMMA)

# Import vector of risk functionals ---------------------------------------
################
l_name<-c("U","Surcote","Hs")
END<-length(l_name)+1
Name_import<-"residuals_MV/RiskF_"
L<-length(l_name)-1
for(z in c(1:L)){
  name_z<-l_name[z]
  Name_import<-paste0(Name_import,l_name[z],"_")
}
Name_import<-paste0(Name_import,l_name[length(l_name)],".csv")
Vect_l_orig<-read.csv(file = Name_import)[,c(2:END)]

# Hidden regular variations ? ---------------------------------------------
#########
pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_min.png"))
Min_Risk_functionals_ALL<-apply(X =  Vect_l_orig,MARGIN = 1,
                                FUN = min)
Vect_k<-c(10:300)
LIMS_Y<-c(0,2)
DIMS_plot<-c(15,15,14,14)
Graphics_estimators_gamma(series = Min_Risk_functionals_ALL,
                          vect_k =Vect_k,
                          Title_graphic =" ",
                          dims_elt_text = DIMS_plot,
                          y_lims = LIMS_Y)
dev.off()
#### Max  
pdf(file = paste0("graphiques_MV/Evol_gamma/Evol_gamma_max.png"))
I<-1
J<-2
Max_Risk_functionals_ALL<-apply(X =  Vect_l_orig,MARGIN = 1,
                                FUN = max)
Graphics_estimators_gamma(series = Max_Risk_functionals_ALL,
                          vect_k =Vect_k,
                          Title_graphic =" ",
                          dims_elt_text = DIMS_plot,
                          y_lims = LIMS_Y)
dev.off()


source("functions_MV/Hidden_RV_functions.R")
### Min/Max gamma parameters (HRV assumptions)
Run_diagnostics_Gamma_G(LIMS_Y = LIMS_Y,
                        dims_elt_text = DIMS_plot,
                        Vectors_HTAIL = Vect_l_orig,
                        I = I,J=J,y_lims = LIMS_Y,
                        Vect_k = Vect_k,NAME_Vars = l_name,
                        q = 0.90)
ROOT_analysis<-paste0(getwd(),'/graphiques_MV/')
ROOT_analysis
# Same thing for vector of risk functions ---------------------------------
png(paste0(ROOT_analysis,"asymp_dependencies/riskF_Dep_",
           l_name[I],"_",l_name[J],".png")
    ,width = 1200, height = 600)
par(mfrow = c(1, 2),  # 1 row, 2 plots
    cex.lab = 1,    # axis label size
    cex.axis = 1,   # axis tick label size
    cex.main = 1,   # main title size
    cex = 1)
POT::chimeas(data = cbind(Vect_l_orig[,I],
                          Vect_l_orig[,J])
             ,which = 1,ask = FALSE,
             cex.lab=1.5)
abline(h=0,col="red")
POT::chimeas(data = cbind(Vect_l_orig[,1],
                          Vect_l_orig[,3])
             ,which = 2,ask = FALSE,
             cex.lab=1.5,
             ylabs = rep(NA,2))
mtext(expression(bar(chi)),cex=1.5,las=2,
      outer=TRUE,line=-20.5)
dev.off()
par(mfrow=c(1,1))

######## Select a subvector if needed.
#################
sub_index<-c(1,2,3)
Vect_l_function<-Vect_l_orig[,sub_index]
Th_graph<-exp(3)
ROOT_analysis<-paste0(getwd(),'/graphiques_MV/')
png(filename = paste0(ROOT_analysis,
                      "/asymp_dependencies/Summary_RiskF_ext_depend.png"),
    width = 1000,height=300)
par(mfrow=c(1,ncol(Vect_l_function)))
l_couples<-t(utils::combn(x = c(1:ncol(Vect_l_function)),m = 2))
for(z in c(1:nrow(l_couples))){
  elt<-l_couples[z,]
  i<-elt[1]
  j<-elt[2]
  name_i<-l_name[i]
  name_j<-l_name[j]
  Sub_vect<-Vect_l_function[,c(i,j)]
  Min_line<-apply(Sub_vect,MARGIN = 1,FUN = min)
  Max_line<-apply(Sub_vect,MARGIN = 1,FUN = max)
  Col_graph<-ifelse(Max_line>Th_graph,no = "black",yes = 
                      ifelse(Min_line>Th_graph,yes = "green",
                             no="orange"))
  plot(Sub_vect[,1],Sub_vect[,2],log="xy",
       xlab=name_i,ylab=name_j,col=Col_graph,cex.lab=1.3)
  abline(h=Th_graph,v=Th_graph,col="red")
}
par(mfrow=c(1,1))
dev.off()

#####
delta_moment<-0.5
value_found<-colMeans(Vect_l_function^(delta_moment))
if(min(Vect_l_function[,1])>1){
  value_std_pareto<--(delta_moment-1)^(-1)
  scale_pareto_d<-(value_found/value_std_pareto)^(1/delta_moment)
  Vect_update<-Vect_l_function/scale_pareto_d
}

if((min(Vect_l_function[,1])>0)&&(min(Vect_l_function[,1])<1)){
  value_std_frechet<-factorial(-delta_moment)
  scale_frechet_d<-(value_found/value_std_frechet)^(1/delta_moment)
  Vect_update<-Vect_l_function/scale_frechet_d
}
summary(Vect_update)
## Lower shape parameter --> HRV is likely. 
LROOT_<-"../../ext_independence/estim_dependence/Penalized_least-squares_estimator/"

source(paste0(LROOT_,"Functions/param_estim.R"))
source(paste0(LROOT_,"Functions/Starting_points.R"))
source(paste0(LROOT_,"Functions/Help_functions.R"))
source(paste0(LROOT_,"Functions/cross_validation.R"))

q_norm<-(1)
W_chosen<-c(1,1)

# Mixture models can only be used if g is a modulus -----------------------
# If we have unequal weights, we end up with g not a modulus. 
Name_riskF<-"Weighted_NQ"
# if(Name_riskF=="sum"){
#   Risk_f<-function(x){return(sum(x))}
# 
# }
# if(Name_riskF=="sum_penalized"){
#   W_chosen<-c(0.7,0.3)
#   Risk_f<-function(x){return(sum_penalized(x = x,vect_w = c(0.7,0.3)))}
# }
if(Name_riskF=="Weighted_NQ"){
  Risk_f<-function(x){
    return(Norm_q(x=x,q=q_norm,
                  weights_nq = c(1,1)))
  }
}
I_ext_indep<-2
Graphics_estimators_gamma(series =Vect_update[,I_ext_indep],
                          vect_k =Vect_k,
                          Title_graphic = " ")
J_ext_indep<-3
Graphics_estimators_gamma(series =Vect_update[,J_ext_indep],
                          vect_k =Vect_k,
                          Title_graphic = " ")
Sub_vect_analysed<-Vect_update[,c(I_ext_indep,J_ext_indep)]
RF_data<-apply(X = Sub_vect_analysed,
               MARGIN = 1,FUN = Risk_f)
TS_param<-0.30
k_metric<-mindist_update(data =RF_data,ts = TS_param,method = "ks")
plot(k_metric$Nb_k,k_metric$value_metric,xlab="Number of exceedances",
     ylab="score_metric",type="l")
Graphics_estimators_gamma(series =RF_data,
                          vect_k =Vect_k,
                          Title_graphic = " ")
k <- 0.35
Nb_exceed<-k*nrow(Vect_l_function)

Threshold<-quantile(RF_data,probs =1-k )
Indexes_extremes<-which(RF_data>Threshold)

# If g is a modulus, we have ----------------------------------------------
GAMMA_G<-1
# if(Name_riskF==c("max"){
#   GAMMA_G<-1
# }else{
#   ML_extRemes<-function_ML_extRemes(data_d=RF_data,typeML="GP",
#                                     NB_years=37,k = Nb_exceed)
#   GAMMA_G<-ML_extRemes$estimator
# }

Inds_extremes<-Sub_vect_analysed[Indexes_extremes,]

### Repartition of points of extreme obs (max)
png(filename = paste0(ROOT_analysis,
                      "/asymp_dependencies/Extreme_obs_RiskF.png"),
    width = 1000,height=300)
Nb_pairs<-factorial(ncol(Sub_vect_analysed))/2*factorial(ncol(Sub_vect_analysed)-1)
par(mfrow=c(1,ncol(Sub_vect_analysed)))
l_couples<-t(utils::combn(x = c(1:ncol(Sub_vect_analysed)),m = 2))
for(z in c(1:nrow(l_couples))){
  elt<-l_couples[z,]
  i<-elt[1]
  j<-elt[2]
  name_i<-l_name[i]
  name_j<-l_name[j]
  Sub_vect<-Inds_extremes[,c(i,j)]
  Min_line<-apply(Sub_vect,MARGIN = 1,FUN = min)
  Max_line<-apply(Sub_vect,MARGIN = 1,FUN = max)
  Col_graph<-ifelse(Max_line>Threshold,no = "black",yes = 
                      ifelse(Min_line>Threshold,yes = "green",
                             no="orange"))
  plot(log(Sub_vect[,1]),log(Sub_vect[,2]),
       xlab=name_i,ylab=name_j,cex.lab=1.3,
       col=Col_graph)
  abline(h=log(Threshold),v=log(Threshold),col="red")
}
par(mfrow=c(1,1))
dev.off()

# Gives an input the measure of -------------------------------------------
# empirical stdfEmp on the ranks. ------------------------------------------------------
N<-nrow(Sub_vect_analysed)
d<-ncol(Sub_vect_analysed)
r_init<-ncol(Vect_update)
lambda_init<-10^(-3)
Model_used<-"SSR_row_log"
if(Model_used=="SSR_row_log"){
  points<-c(0,1/4,1/3,1/2,3/4,1)
}else{
  points <- c(0 , 1/6,  1/8 ,  1/4 ,  1/3 , 1/2 , 2/3 , 3/4 , 1)
}
### Non zero is an important parameter
NONZERO<-d-1
Grid_points <- tailDepFun::selectGrid(points, d = d,
                                      nonzero = NONZERO)
R <- apply(Sub_vect_analysed , 2 , rank)
start <- starting_point(Sub_vect_analysed, r_init)
start <- c(t(start))
p <- 0.3
q <- nrow(Grid_points)
w_total <- sapply(1:q, function(m)
  tailDepFun::stdfEmp(R, k * N, Grid_points[m, ]))
# Use functions from mixture models ---------------------------------------
####################
dyn.load("main.dll")
for(Number_mod in c(1:5)){
  start <- starting_point(Vect_update, Number_mod)
  start <- c(t(start))
  Optimization<-param_estim(d =d, r = Number_mod, grid = Grid_points , 
                            lambda=lambda_init , 
                                  num_col = NULL ,
                                  start = start , type = Model_used,
                                  p = p ,  w = w_total )
}
R_chosen<-2
start_Rchosen <- starting_point(Vect_update, R_chosen)
start_Rchosen <- c(t(start_Rchosen))
Params_first_estim<-param_estim(d =d, r = R_chosen, grid = Grid_points , 
            lambda=lambda_init , 
            num_col = NULL ,
            start = start_Rchosen, type = Model_used,
            p = p ,  w = w_total )
Params_first_estim$pls_matrix
Theta_found<-Params_first_estim$pls_dep
Theta_found
### The maximum of the second variable (surge)
### is only modeled by the fourth model,
### which only intervenes with this variable. 
### While the two other ones are associated with 
### two models.
require(tailDepFun)
K<-5
w_data <- W_calculus(k = k * N , num_class = K , X = Sub_vect_analysed, 
                grid = Grid_points, q = q)
#Vect_lambda<-c(1e-04,1e-03,1e-02,1e-01)
Vect_lambda <- seq(10^(-6), 1 , length.out = 50)
L_Lambda<-list()
Z<-1

# To improve, does not depend on lambda -----------------------------------

for (lam in Vect_lambda){
  CV_run<-cross_validation(d , R_chosen,  grid = Grid_points,
                           lambda = lam , num_col = NULL, start = start_Rchosen ,
                           type = Model_used,
                           p = p, w_data , num_class=K)
  Optimization_lam<-param_estim(d =d, r =R_chosen, grid = Grid_points , 
                                lambda=lam, 
                            num_col = NULL ,
                            start = start_Rchosen, type = Model_used,
                            p = p ,  w = w_total )
  L_Lambda[[Z]]<-list("cv"=CV_run,"result_theta_A"=Optimization_lam)
  Z<-Z+1
}
Result_CV<-sapply(L_Lambda,function(x){return(x$cv)})
plot(Vect_lambda,Result_CV,
     main="CV score according to lambda value")
Index_lambda_opt<-which.min(Result_CV)
Params_opt<-L_Lambda[[Index_lambda_opt]]

# Simulation of the T-generator -------------------------------------------
###########
ROOT_import_simul<-"../../ext_independence/work_mixture/MGPD-simulation/"

source(paste0(ROOT_import_simul,"Helper_functions.R"))
## Convergence of angle measure of the g Pareto process 
Aopt<-Params_opt$result_theta_A$pls_matrix
Aopt
Theta_opt<-Params_opt$result_theta_A$pls_dep
Theta_opt
if(Model_used=="SSR_row_HR"){
  diag(Theta_opt)<-rep(0,ncol(Theta_opt))
  Sigma<-list()
  Theta_opt<-Params_opt$result_theta_A$pls_dep
  diag(Theta_opt)<-rep(0,ncol(Theta_opt))
  for(k in c(1:R_chosen)){
    Sigma_tilde<-GAMMA_to_SigmaTilde(Mat_GAMMA = Theta_opt,k_const=k)
    Sigma[[k]]<-Sigma_tilde
  }
}


Summary_A_coefficients<-Aopt
rownames(Summary_A_coefficients)<-names(Sub_vect_analysed)
Summary_A_coefficients<-round(Summary_A_coefficients,2)
Summary_A_coefficients
write.csv(x=Summary_A_coefficients,
    file = paste0(getwd(),"/residuals_MV/Latent_repartition.csv"))
# Var(Z_k)=0, Cov(Z_k,Z_j)=0 as we have a constant
# comes from the definition of Hussler-Reiss (See Graphical models by Engelke)

# Test for simulations.

M_simul<-1000
Model_used
# For simulations, change the shape parameter if the 
# risk function is not the maximum. 
Theta_opt
Aopt
if(Model_used=="SSR_row_log"){
  source(paste0(ROOT_import_simul,"mgpd_simulation_mixture_logistic.R"))
  Result_sim<-t(replicate(n = M_simul,mgpd_simulation_mixture_logistic_modif(d =d,
                                                   r = R_chosen, alpha =rep(Theta_opt,
                                                                            R_chosen),
                                                   A = Aopt,RiskF = max,
                                                  Shape_g= GAMMA_G)))
}else{
  require(mvtnorm)
  source(paste0(ROOT_import_simul,"mgpd_simulation_mixture_HR.R"))
  Result_sim<-t(replicate(n = M_simul,mgpd_simulation_mixture_HR_modif(d =d,
                                                                 r = R_chosen,
                                                                 A = Aopt,
                                                                 Sigma = Sigma,
                                                                 Shape_g= GAMMA_G)))
}
Risk_f
vect_extremes<-Sub_vect_analysed[Indexes_extremes,]
Result_sim_MxF<-exp(Result_sim)*Threshold
summary(vect_extremes)
summary(Result_sim_MxF)
png(filename = paste0(ROOT_analysis,
                      "/asymp_dependencies/Compar_SimulMixt_RiskF.png"),
    width = 1000,height=300)
Nb_pairs<-factorial(ncol(Result_sim))/2*factorial(ncol(Result_sim)-2)
par(mfrow=c(1,Nb_pairs))
l_couples

for(elt in l_couples){
  i<-elt[1]
  j<-elt[2]
  name_i<-l_name[i]
  name_j<-l_name[j]
  print(c(i,j))
  plot(vect_extremes[,i],vect_extremes[,j],
       log="xy",xlab=name_i,ylab=name_j)
  points(Result_sim_MxF[,i],Result_sim_MxF[,j],col="red",
         pch=20)
}
par(mfrow=c(1,1))
dev.off()
Sigma
MAX_law<-apply(X = Result_sim,MARGIN = 1,FUN = min)
goftest::ad.test(x = MAX_law,null = "pexp",rate=1)
### Does follow an exponential law.
# Exp to Pto
EXP_to_PTO<-(exp(Result_sim))
Graphics_estimators_gamma(series = EXP_to_PTO[,3],
                          vect_k = Vect_k,
                          Title_graphic = " ")
PTO_question<-apply(X = EXP_to_PTO,MARGIN = 1,
                    FUN = max)
Graphics_estimators_gamma(series = PTO_question,
                          vect_k = Vect_k,
                          Title_graphic = " ")
goftest::ad.test(x =PTO_question,null = extRemes::"pevd",
                 shape=1,threshold=1,scale=1,type="GP")
# Frech to Pto
goftest::ad.test(x=RF_data[Indexes_extremes],null = extRemes::"pevd",
                 shape=1,threshold=as.numeric(Threshold),
                 scale=as.numeric(Threshold),type="GP")
### Plot simulations obtained
png(filename = paste0(ROOT_analysis,
                      "/asymp_dependencies/SimulMixt_RiskF.png"),
    width = 1000,height=300)
par(mfrow=c(1,ncol(Result_sim)))
for(elt in l_couples){
  i<-elt[1]
  j<-elt[2]
  name_i<-l_name[i]
  name_j<-l_name[j]
  Sub_vect<-Result_sim[,c(i,j)]
  Min_line<-apply(Sub_vect,MARGIN = 1,FUN = min)
  Max_line<-apply(Sub_vect,MARGIN = 1,FUN = max)
  Col_graph<-ifelse(Max_line>0,no = "black",yes = 
                      ifelse(Min_line>0,yes = "green",
                             no="orange"))
  plot(Sub_vect[,1],Sub_vect[,2],
       xlab=name_i,ylab=name_j,cex.lab=1.3,
       col=Col_graph)
  abline(h=0,v=0,col="red")
}
par(mfrow=c(1,1))
dev.off()

# Non param ---------------------------------------------------------------
#### ---------------------------------------------------------------------

EXPS<--log(1-UNIFS)
summary(EXPS)
L_exp<-apply(X = EXPS,MARGIN = 1,FUN = sum)
Seuil_lg<-quantile(L_exp,0.75)
Inds_extremes<-which(L_exp>Seuil_lg)
# plot(density(EXPS[Inds_extremes,1]-EXPS[Inds_extremes,2]))
# Result_obtained<-Non_param_LEGRAND_RiskF(Data_scale_exp = EXPS[Inds_extremes,],
#                         nb_simul = 1000,MOD_EGPD_RiskF =MOD_EGPD_RISKF,
#                         Threshold_EXP = Seuil_lg)
UNIFS
AD_test<-goftest::ad.test(x = UNIFS,null = "punif")$p.value
AD_test
KS_test<-ks.test(x = UNIFS,"punif")$p.value
KS_test
CVM_test<-goftest::cvm.test(UNIFS,"punif")$p.value
CVM_test

# Param -------------------------------------------------------------------
PARETO<--1/(log(UNIFS))
plot(PARETO,log="xy")

#Use scoring rule to choose the correct params ---------------------------
############### ----------------------------------------------------------
L_g<-apply(X = Vect_l_function,MARGIN = 1,
           FUN = function(x){return(sum(x))})
# L_g<-apply(X = PARETO,MARGIN = 1,FUN = sum)
Seuil_lg<-quantile(L_g,0.75)
Seuil_lg
Inds_exts<-which(L_g>Seuil_lg)
plot(PARETO[Inds_exts,],log="xy")
abline(h=Seuil_lg,v=Seuil_lg,col="red")
POT::tcplot(L_g[Inds_exts],ask = FALSE)

# Weight function ---------------------------------------------------------
# weightFun <- function(x, u){
#   SUM<-sum(x / u) - 1
#   return(x * (1 - exp(-(SUM))))
# }
# 
# #Define partial derivative of weighting function
# dWeightFun <- function(x, u){
#   (1 - exp(-(sum(x / u) - 1))) + (x / u) * exp( - (sum(x / u) - 1))
# }

ALPHA<-0.5
Exts_list <- split(PARETO[Inds_exts,], 
                   seq(nrow(PARETO[Inds_exts,])))
source("functions_MV/function_estim_simul_gPareto.R")
Objective_function(weightFun = weightFun_sum,
                   dWeightFun = DweightFun_sum,Extreme_inds  =Exts_list,
                   u=Seuil_lg,theta =0.5,model_used = "log" )
Opt_reduc<-optim(par =0.5,
      fn = Objective_function,
      method = "L-BFGS-B",weightFun = weightFun_sum ,
      dWeightFun = DweightFun_sum,Extreme_inds  =Vect_l_function[Inds_exts,],
      u=Seuil_lg,lower = 0,upper = 1,
      model_used="log")
Opt_reduc$par
Result_alpha<-function(ALPHA_test){
  return(Objective_function(ALPHA = ALPHA_test,weightFun = weightFun ,
                            dWeightFun = dWeightFun,Extreme_inds  =Exts_list,
                            u=Seuil_lg))
}

# Simulation of the couple of l for general g -------------------------------------------

Test<-Rejection_sampling_ParetoP(params_model = 1/ALPHA_opt_scoring,
                           name_model = "log",
                           riskF_defined = sum, 
                           Threshold = Seuil_lg,
                           d = 2,sum_used = TRUE,
                           M = 1000)

# Compare simulations and obsevations ----------------------------------------------------
##### --------------------------------------------------------------------
plot(PARETO[Inds_exts,],log="xy",col="blue")
abline(h=Seuil_lg,v=Seuil_lg,col="red")
points(Test,col="orange")

summary(PARETO[Inds_exts,])
summary(Test)
Value_found<-apply(X = Test,MARGIN = 1,FUN = sum)
summary(Value_found)

goftest::ad.test(x = Value_found,null = extRemes::"pevd",
                 shape=1,scale=Seuil_lg,threshold=Seuil_lg)

# Normalisation constants of the maximum -------------------------------------------------
#######-----------------------------------

Nb_years<-37
NPobs<-nrow(Vect_l_function)/Nb_years
U<-20
V<-U
Result<-1-as.numeric((Vect_l_function[,1]/NPobs<U)&(Vect_l_function[,2]/NPobs<V))
Lambda<-sum(Result)/Nb_years
Lambda
Dates<-read.csv(file = "../ss_tend/HIVER/dates_prises.csv")$x
Translated_dates<-as.POSIXct(Dates, format="%d/%m/%Y")
Decimal_date<-lubridate::decimal_date(Translated_dates)
#from diplot
date <- floor(Decimal_date)
tim.rec <- range(Decimal_date)
nb.occ <- NULL
for (year in tim.rec[1]:tim.rec[2]) nb.occ <- c(nb.occ, 
                                                sum(Result & (date == year)))
nb.occ
Variance_<-var(nb.occ)
DIndex<-Lambda/Variance_
conf<-0.95
conf_sup <- qchisq(1 - (1 - conf)/2, Nb_years - 1)/(Nb_years - 1)
conf_inf <- qchisq((1 - conf)/2, Nb_years- 1)/(Nb_years - 1)
Vect_scale<-log(1+Vect_l_function[,1])

# Approaching the marginals -----------------------------------------------
Kappa<-1
summary(Vect_l_function)
Vect_puissance<-log(Vect_l_function+1)
Threshold<-log(20+1)
MODEV<-mev::gp.fit(xdat = Vect_puissance[,1],threshold = Threshold)
Ajust_EGP_1<-mev::fit.extgp(Vect_puissance[,1],
                        model = 1,method = "mle",
                        init = c(Kappa,MODEV$est),
                        R = 20)
Fit_EGP_1<-Ajust_EGP_1$fit$mle
MODEV2<-mev::gp.fit(xdat = Vect_puissance[,2],threshold = Threshold)
Ajust_EGP_2<-mev::fit.extgp(Vect_puissance[,2],
                            model = 1,method = "pwm",
                            init = c(Kappa,MODEV2$est),
                            R = 20)
Fit_EGP_2<-Ajust_EGP_2$fit$pwm

# Conversion to Unif ------------------------------------------------------
Sub<-which(Vect_puissance[,1]>0)
Unif_1<-sapply(Vect_puissance[Sub,1],function(x){
  return(mev::pextgp(q =x ,kappa =Fit_EGP_1[["kappa"]],
                    xi = 0,sigma = Fit_EGP_1[["sigma"]],
                    type = 1))})
plot(density(Unif_1))
summary(Unif_1)
goftest::ad.test(x = Unif_1)

Fill<-matrix(NA,ncol=ncol(Vect_input),nrow=3)
for(j in c(1:ncol(Vect_input))){
  Result_std[,j]<-Vect_input[,j]/scale_frechet_d[j]
  pval<-goftest::ad.test(Vect_input[,j],null = extRemes::"pevd",
                         scale=scale_frechet_d[j],shape=1,
                         loc=scale_frechet_d[j],type="GEV")$p.value
  pval2<-ks.test(Vect_input[,j],extRemes::"pevd",
                 scale=scale_frechet_d[j],shape=1,
                 loc=scale_frechet_d[j],type="GEV")$p.value
  pval3<-goftest::cvm.test(Vect_input[,j],extRemes::"pevd",
                           scale=scale_frechet_d[j],shape=1,
                           loc=scale_frechet_d[j],type="GEV")$p.value
  Fill[,j]<-c(pval,pval2,pval3)
}

typemod<-"log"
SYMbool<-FALSE
AIC_ext_d_model<-function(data,method_opt,model,SYMBOOL){
  Mod_biv<-evd::fbvevd(x = data,
                       method=method_opt,model=model,
                       sym=SYMBOOL,std.err = FALSE)
  return(AIC(Mod_biv))
}
Families<-c("log","alog","hr","bilog","ct",
            "neglog","amix")
L_AIC<-sapply(Families,AIC_ext_d_model,
       method_opt="BFGS",
       SYMBOOL=FALSE,data=Result_std)
L_AIC
Familyopt<-Families[which.min(L_AIC)]
Familyopt
##### (hr) model seems to be slightly better. 
Mod_biv_optAIC<-evd::fbvevd(x = Result_std,
            method="BFGS",model=Familyopt,
            sym=FALSE,std.err = FALSE)
plot(Mod_biv_optAIC,ask = FALSE)

# Test statistic param vs non param----------------------------------------------------------
delta<-0.05
w<-seq.int(0, 1,by=delta)
Anonparam<-evd::abvnonpar(w,data =Result_std)
plot(w,Anonparam,type="l")
Aparam<-evd::abvevd(w,dep = Mod_biv_optAIC$estimate[["dep"]],
                    asy = Mod_biv_optAIC$estimate[["asy1"]] ,
                    beta = Mod_biv_optAIC$estimate[["beta"]],
                    alpha = Mod_biv_optAIC$estimate[["alpha"]],
                    model = Mod_biv_optAIC$model)
Sn_Distance<-sum(nrow(Result_std)*(Aparam-Anonparam)**(2)*delta)
Sn_Distance

MAR1<-c(Mod_biv_optAIC$estimate[["loc1"]],
        Mod_biv_optAIC$estimate[["scale1"]],
        Mod_biv_optAIC$estimate[["shape1"]])
MAR2<-c(Mod_biv_optAIC$estimate[["loc2"]],
        Mod_biv_optAIC$estimate[["scale2"]],
        Mod_biv_optAIC$estimate[["shape2"]])
# Chi plot vs model -------------------------------------------------------
Chi_Value_craft<-2+log(evd::pbvevd(q=c(1,1),
            model = Mod_biv_optAIC$model,
            dep=Mod_biv_optAIC$estimate[["dep"]],
            mar1 = MAR1,
            mar2 = MAR2,
            alpha=Mod_biv_optAIC$estimate[["alpha"]],
            beta = Mod_biv_optAIC$estimate[["beta"]],
            asy = Mod_biv_optAIC$estimate[["asy"]]))
Chi_Value<-Mod_biv_optAIC$dep.summary
par(mfrow=c(1,2))
POT::chimeas(Result_std,which=1)
abline(h=Chi_Value,col="red")
abline(h=0,col="black")
POT::chimeas(Result_std,which=2)
par(mfrow=c(1,1))

source("residus_MV/Stat_test_MV.R")
Boot_rep<-200
Distrib_Sn_pb<-Stat_Test_Measure(Boot_rep = Boot_rep,Mod_BIV =Mod_biv_optAIC,
                                 Size_data = nrow(Result_std),
                                 SYMbool = SYMbool,
                                 w = w,MAR1=MAR1,MAR2=MAR2)
plot(density(Distrib_Sn_pb),
     main="Distribution of parametric bootstrap statistic")
abline(v=Sn_Distance,col="red")
summary(Distrib_Sn_pb)

Lambda_model<--log(evd::pbvevd(q=c(U,V),
                model = Mod_biv_optAIC$model,
                dep=Mod_biv_optAIC$estimate[["dep"]],
                mar1 = MAR1,
                mar2 = MAR2,
                alpha=Mod_biv_optAIC$estimate[["alpha"]],
                beta = Mod_biv_optAIC$estimate[["beta"]],
                asy = Mod_biv_optAIC$estimate[["asy"]]))
Lambda_model

# Data point process ------------------------------------------------------
plot(cumsum(Result),type="l",ylab="Immerged chain",xlab="Time")

# Model point process -----------------------------------------------------
Nb_jumps<-max(cumsum(Result))
Time_jump<-Generator_PP_wake_up(Number_realisations = Nb_jumps,
                                lambda = Lambda_model)
Time_maxunit<-Time_jump*NPobs
points(round(Time_maxunit),c(1:Nb_jumps),col="red",
       type = "b")
