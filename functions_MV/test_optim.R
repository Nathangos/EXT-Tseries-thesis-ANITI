rm(list=ls())
set.seed(133)
require(mev)
require(reshape2)
require(ggplot2)
source("../fonctions/fonctions.R")
source("functions_MV/function_estim_simul_gPareto.R")
source("functions_MV/function_simul_MV.R")

# TEST difference models --------------------------------------------------
# Risk_f<-function(x){return(sum(x))}
#Name_riskF<-"NQ"
W_chosen<-c(0.7,0.3)
### Max/q<-40
### Min/q<--40
q<--2
Name_riskF<-"Weighted_NQ"
if(Name_riskF=="sum"){
  Risk_f<-function(x){
    return(sum(x))}
  wFun<-weightFun_sum
  dwFun<-DweightFun_sum

}
if(Name_riskF=="sum_penalized"){
  Risk_f<-function(x){
  return(sum_penalized(x = x,vect_w = W_chosen))}
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
  Risk_f<-function(x){
    return(Norm_q(x=x,q=q,
                  weights_nq = W_chosen))
  }
  wFun<-function(x,u){
    return(weightFun_nq(x = x,u = u,
                        q = q,
                        weights_nq = W_chosen))
  }
  dwFun<-function(x,u){
    return(DweightFun_nq(x = x,u = u,
                         q = q,weights_nq = W_chosen))
  }
}

d<-2
Define_Gamma_from_lambda<-function(lambda,d){
  dalt<-d-1
  ELT_base<-rep(c(0,rep(lambda,d)),dalt)
  mat_gamma<-matrix(ELT_base,nrow = d,ncol=d,
                     byrow = FALSE)
  return(mat_gamma)
}
model_used<-"log"
if(model_used=="log"){
  THETA<-0.8
}
if(model_used=="hr"){
  THETA<-Define_Gamma_from_lambda(lambda = 0.5,d = d)
}
THETA
if(model_used=="log"){
  MEV_FRECHET<-rmev(n=5000, d=d, param=1/ THETA, model='log', alg='ef')
  INIT_THETA<-0.3
}else{
  MEV_FRECHET<-rmev(n=5000, d=d, sigma = THETA, model='hr')
  INIT_THETA<-Define_Gamma_from_lambda(lambda = 0.3,d = d)
  INIT_THETA<-INIT_THETA[t(utils::combn(x = c(1:d),m = 2))]
}
Sub_theta<-INIT_THETA

List_Params_RF<-list("parametric"=TRUE,"name_RF"=Name_riskF,
                     "function"=Risk_f,"name_model"=model_used,
                     "Scale_Frechet"=c(1,1),"weightFunction"=wFun,
                     "dWeightFunction"=dwFun,"init_opt_param"=Sub_theta)
List_Params_RF
summary(apply(MEV_FRECHET,MARGIN = 1,FUN = min))
Sum_found<-apply(X = MEV_FRECHET,MARGIN = 1,FUN = Risk_f)

summary(Sum_found)
TS_param<-0.50
k_metric<-mindist_update(data = Sum_found,ts = TS_param,
                         method = "ks")
kZero_found<-which.min(k_metric$value_metric)
plot(k_metric$Nb_k,k_metric$value_metric,xlab="Number of exceedances",
     ylab="score_metric",type="l")
abline(v=kZero_found,col="red")
# Graphics_estimators_gamma(series = Sum_found,
#                           vect_k = c(10:500),
#                           Title_graphic = " ",
#                           NB_years = NULL)
Kfound<-200
Q<-Kfound/length(Sum_found)
Th<-quantile(x = Sum_found,probs = 1-Q)
Risk_f_q<-function(x,Q_var){
  return(Norm_q(x=x,q=Q_var,
                weights_nq = W_chosen))
}
VEct_q<-seq.int(from = -10,to = 10,by = 1)
Vect_diff_Q<-sapply(X = VEct_q,FUN=function(q){
  Vect_for_oneQ<-apply(X = MEV_FRECHET,MARGIN = 1,FUN=function(x,Q){
    return(Risk_f_q(x =x,Q_var=q))
  },Q=q)
  return(Vect_for_oneQ)
})
V_diffQ<-as.data.frame(Vect_diff_Q)
colnames(V_diffQ)<-as.character(V_diffQ)
THETA_found<-Estim_param_RF_homogeneous(Params_risk_Function = List_Params_RF,
                           Seuil_lprime =Th ,
                           Vect_l_function =MEV_FRECHET,
                           d=d)
THETA_found
Gap<-abs(THETA_found-THETA)
print(paste0("The absolute difference between the estimator and the true value is ",Gap))
Inds<-which(Sum_found>Th)
Exts_mev<-MEV_FRECHET[Inds,]
if(sum(W_chosen)==1){
  png(filename = paste0("functions_MV/Weighted_results_Lp=",q,".png"))
  plot(MEV_FRECHET,log="xy",xlab="First dimension",ylab="Second dimension",
       cex.lab=1.5)
  abline(h=Th,v=Th,lty="dashed",col="blue")
}else{
  png(filename = paste0("functions_MV/results_Lp=",q,".png"))
  plot(MEV_FRECHET,log="xy",xlab="First dimension",ylab="Second dimension",
       cex.lab=1.5)
  abline(h=Th,v=Th,lty="dashed",col="blue")
}

## Change the threshold used.
# J_plot<-2
# X_private<-Exts_mev[,-J_plot]
# seq_x<-seq.int(from = 0.001,to =20,length.out =100)
# seq_xj<-sapply(seq_x,FUN=Shape_boundary,
#                j = J_plot,ul = Th,RiskF = Risk_f)
# lines(seq_x,seq_xj,col="blue")
M<-3000
#### Directly gives the lambda square coefficients for rgpapr != mev, gives the sigma...
Simuls_gp<-Simulation_gParetoP(Params_risk_Function=List_Params_RF,
                    Threshold= Th,d=ncol(MEV_FRECHET),
                    M=M,theta_opt=THETA_found)
THETA_found
stats_obtained<-apply(Simuls_gp,MARGIN = 1,FUN = Risk_f)
points(Simuls_gp,col="red",
       cex=0.5)
dev.off()

if(List_Params_RF[["name_RF"]]!="max"){
  Evol_cost<-function(theta_candidate){
    Objective_function(theta = theta_candidate,
                       weightFun = List_Params_RF[["weightFunction"]],
                       dWeightFun = List_Params_RF[["dWeightFunction"]],
                       Extreme_inds =as.list(as.data.frame(t(Exts_mev))) ,
                       u = Th,
                       model_used = List_Params_RF[["name_model"]])
    
  }
  if(model_used=="log"){
    vect_candidates<-seq.int(from = 0.1,to = 0.99,length.out = 20 )
  }
  if(model_used=="hr"){
    vect_lambda<-seq.int(from = 0.1,to = 0.99,length.out = 20)
    vect_candidates<-lapply(vect_lambda,Define_Gamma_from_lambda,d = d)
  }
  Cost_evol<-sapply(vect_candidates,
                    Evol_cost)
  plot(vect_candidates,Cost_evol)
}

# Test mixture model ------------------------------------------------------
####
Name_riskF<-"sum_penalized"
#Name_riskF<-"sum"
if(Name_riskF=="sum"){
  Risk_f<-function(x){return(sum(x))}
  
}
if(Name_riskF=="sum_penalized"){
  W_chosen<-c(0.7,0.3)
  Risk_f<-function(x){return(sum_penalized(x = x,vect_w = c(0.7,0.3)))}
}
NB<-100
Th<-5
Result_test<-t(Mixture_d(prop_bernoulli = c(0.5,0.5),
          nb_simulations =NB,d = 2,RiskFunction=max))
summary(Result_test)
Result_T<-Result_test-apply(X = Result_test,MARGIN = 1,
                  FUN = Risk_f)
Result_T<-Result_T+rexp(n = NB,rate = 1)+Th
summary(Result_T)
### End file of code test.
#########################
