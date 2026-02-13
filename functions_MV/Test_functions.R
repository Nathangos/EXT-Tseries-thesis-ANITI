rm(list=ls())
set.seed(133)
require(mev)
### Source des fonctions. ####
#Functions from univariate case. 
source("../fonctions/fonctions_Pareto.R")
source("../fonctions/fonctions.r")
source("../fonctions/fonctions_r_pareto_EVOLTAU.r")
source("functions_MV/Hidden_RV_functions.R")

Simul_mev<-mev::rmev(n = 1000,d = 2,param = 1,model="log")
Name_vars<-c("Mi",
             "Umi")
Vect_k<-c(10:400)
LIMS_Y<-c(0,2)
DIMS_plot<-c(15,15,14,14)
Xi_<-apply(X = Simul_mev,MARGIN = 1,
           FUN = sum)
Eta_<-Simul_mev[,1]/Xi_
## Theta_1
Hillish_values_1<-sapply(X = Vect_k,FUN = Hillish_stat,
                         Eta=Eta_,Xi=Xi_)
Hillish_values_2<-sapply(X = Vect_k,FUN = Hillish_stat,
                         Eta=(-1)*Eta_,Xi=Xi_)
## Theta_2
Xi_2<-Min_Risk_functionals[-Indexes_theta_1]
Sub_couple<-couple_UV[-Indexes_theta_1,]
Eta_2<-Sub_couple[,1]/Sub_couple[,2]
Hillish_values_minus_1<-sapply(X = Vect_k,FUN = Hillish_stat,
                               Eta=Eta_2,Xi=Xi_2)
Hillish_values_minus_2<-sapply(X = Vect_k,FUN = Hillish_stat,
                               Eta=(-1)*Eta_2,Xi=Xi_2)
par(mfrow=c(2,2))
plot(Vect_k,Hillish_values_1,type="l",
     xlab="Number of exceedances",
     ylab="Hillish value",cex.lab=1.5)
plot(Vect_k,Hillish_values_2,type="l",
     xlab="Number of exceedances",
     ylab="Hillish value",cex.lab=1.5)
plot(Vect_k,Hillish_values_minus_1,type="l",
     xlab="Number of exceedances",
     ylab="Hillish value",cex.lab=1.5)
plot(Vect_k,Hillish_values_minus_2,type="l",
     xlab="Number of exceedances",
     ylab="Hillish value",cex.lab=1.5)
par(mfrow=c(1,1))
dev.off()

# Run_diagnostics_Gamma_G(LIMS_Y =LIMS_Y,
#                         dims_elt_text = DIMS_plot,
#                         Vectors_HTAIL = Simul_mev,
#                         I = 1,J = 2,NAME_Vars = Name_vars,
#                         Vect_k = Vect_k,q=0.8)

source("functions_MV/function_estim_simul_gPareto.R")
source("functions_MV/function_simul_MV.R")
### with param model
TS_test<-0.3
Simul_test<-mev::rmev(n = 10000,d = 2,param = 0.2,model = "log")
POT::chimeas(Simul_test,which = 1)
Mstat<-apply(X = Simul_test,MARGIN = 1,FUN = max)
k_metric<-mindist_update(data =Mstat,
                         ts = TS_test,method = "mad")
plot(k_metric$Nb_k,k_metric$value_metric,xlab="Number of exceedances",
     ylab="score_metric",type="l")
kZero_found<-which.min(k_metric$value_metric)
Qtest<-1-(kZero_found/nrow(Simul_test))
Qtest<-0.90
##Exts_inds
####
Nu<-quantile(Mstat,Qtest)
Inds_exts<-which(Mstat>Nu)
Exts_mev<-Simul_test[Inds_exts,]
POT::chimeas(Exts_mev,which=1)
POT::chimeas(Exts_mev,which=2)

UNIF<-exp(-Simul_test^(-1))
EXP<--log(1-UNIF)
goftest::ad.test(x = EXP[,2],null = "pexp")
Output<-Simul_RF_max_NP(Params_risk_Function = List_Params_RF,
                        M=2000,QTH = Qtest,
                        Vect_l_function = EXP)
Lg<-apply(X =EXP,MARGIN = 1,FUN = max)
Threshold_exp<-quantile(Lg,Qtest)
Output_drift<-pexp(Output+Threshold_exp)
par(mfrow=c(1,2))
POT::chimeas(EXP[Inds_exts,],which=1)
POT::chimeas(Output_drift,which = 1)
par(mfrow=c(1,1))
plot(UNIF[Inds_exts,])
points(Output_drift,col="red")

### test independent gumbel
th<-5
n<-1000
sig1<-3
sig2<-5
U<--log(-log(runif(n = n)))*sig1
V<--log(-log(runif(n = n)))*sig2
cor.test(x = U,y = V)
Tvect<-cbind(U,V)
E<-rexp(n)+th
Sim_MGP<-Tvect+-apply(X = Tvect,MARGIN = 1,FUN = max)+E
plot(Sim_MGP)
cor.test(x = Sim_MGP[,1],y = Sim_MGP[,2])
POT::chimeas(Sim_MGP,which = 1)
### Do find asymptotic dependence

### Normal vs 2 component Gaussian mixture

d<-2
V_dummy<-mvtnorm::rmvnorm(n = 1000,mean = rep(0,d))
TwoGM_versus_normal(vect_Id = V_dummy,
                    alpha = 0.05)
### 
epsi<-0.2
prob_classes<-c(1-epsi,epsi)
nu<-rep(0.3,d)
Test_1sim<-function(prob_classes,nu,n){
  Rep_Normal<-mvtnorm::rmvnorm(n = n,mean = rep(0,d))
  Result_test<-TwoGM_versus_normal(vect_Id = Rep_Normal,alpha = 0.05)
  return(Result_test)
}
Test_procedure<-function(N_test,prob_classes,nu,n){
  vect_rep_summary<-replicate(n = N_test,Test_1sim(prob_classes = prob_classes,
                                                   nu = nu,n = n))
  return(vect_rep_summary)
}
Result_pval<-Test_procedure(N_test = 1000,prob_classes = prob_classes,
                            nu = nu,n = 1000)
Mat_R<-colMeans(t(apply(X = Result_pval,MARGIN = 2,FUN = unlist)))


