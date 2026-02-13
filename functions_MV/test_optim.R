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
q<-20
Name_riskF<-"max"
if(Name_riskF=="max"){
  Risk_f<-function(x){return(max(x))}
  wFun<-NA
  dwFun<-NA
}
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
# # ,
# weights_nq = W_chosen
#,weights_nq = W_chosen
#,weights_nq = W_chosen
if(Name_riskF=="Weighted_NQ"){
  Risk_f<-function(x){
    return(Norm_q(x=x,q=q,weights_nq = c(0.3,0.7)))
  }
  wFun<-function(x,u){
    return(weightFun_nq(x = x,u = u,
                        q = q,weights_nq = c(0.3,0.7)))
  }
  dwFun<-function(x,u){
    return(DweightFun_nq(x = x,u = u,
                         q = q,weights_nq = c(0.3,0.7)))
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
model_used<-"bilog"
if(model_used=="bilog"){
  THETA<-c(0.50,0.55)
}
if(model_used=="log"){
  THETA<-0.80
}
if(model_used=="hr"){
  THETA<-Define_Gamma_from_lambda(lambda = 0.5,d = d)
}
if(model_used=="bilog"){
  print("here")
  MEV_FRECHET<-rmev(n=5000, d=d, param = THETA, 
                    model=model_used, alg='ef')
  INIT_THETA<-c(0.3,0.7)
}
print(THETA)
if(model_used=="log"){
  MEV_FRECHET<-rmev(n=100, d=d, param=1/THETA, model='log', alg='ef')
  INIT_THETA<-0.3
}
summary(MEV_FRECHET)
goftest::ad.test(x =exp(-(MEV_FRECHET[,1])^(-1)),
                 null="punif")
# else{
#   MEV_FRECHET<-rmev(n=5000, d=d, sigma = THETA, model='hr')
#   INIT_THETA<-Define_Gamma_from_lambda(lambda = 0.3,d = d)
#   INIT_THETA<-INIT_THETA[t(utils::combn(x = c(1:d),m = 2))]
# }
RF<-rowSums(MEV_FRECHET^2)^(1/2)
summary(RF)
Graphics_estimators_gamma(series =RF,
                          vect_k = c(70:200),
                          Title_graphic = "test RF of Frechet distrib")
goftest::ad.test(RF,null = extRemes::"pevd",
                 type="GEV",shape=1,scale=2,loc=2)
TEST<-exp(-(MEV_FRECHET)^(-1))
TEST<--log(1-TEST)
MAX<-apply(X = TEST,MARGIN = 1,FUN = max)
SEUIL_EXP<-quantile(MAX,0.98)
INDS_exts<-which(MAX>SEUIL_EXP)
goftest::ad.test(MAX[INDS_exts]-SEUIL_EXP,null = "pexp")
SUB<-TEST[INDS_exts,]
Simuls_<-Non_param_LEGRAND_RiskF(Data_scale_exp = SUB,
                        nb_simul = 300,Threshold_EXP = SEUIL_EXP)
GAP<-SUB
Delta_i<-as.numeric(t(diff(t(GAP))))*(-1)
head(Delta_i)
head(GAP)
nb_simul<-300
Delta_tilde<-sample(Delta_i,size = nb_simul,
                    replace = TRUE)
IND_plus<-which(Delta_tilde>0)
D<-2
Delta_mat<-matrix(0,ncol = D,
                  nrow=nb_simul)
Delta_mat[-IND_plus,1]<-Delta_tilde[-IND_plus]
Delta_mat[IND_plus,2]<--Delta_tilde[IND_plus]
Simul_Z<-rexp(n=nb_simul)
Simul_RiskF_value<-Simul_Z+Delta_mat+SEUIL_EXP
plot(SUB)
points(Simul_RiskF_value,col="red")

List_Params_RF<-list("parametric"=TRUE,"name_RF"=Name_riskF,
                     "function"=Risk_f,"name_model"=model_used,
                     "Scale_Frechet"=c(1,1),"weightFunction"=wFun,
                     "dWeightFunction"=dwFun,"init_opt_param"=INIT_THETA)
List_Params_RF
summary(apply(MEV_FRECHET,MARGIN = 1,FUN = min))
Sum_found<-apply(X = MEV_FRECHET,MARGIN = 1,FUN = Risk_f)

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
Kfound<-150
Q<-1-(Kfound/length(Sum_found))
Th<-quantile(x = Sum_found,probs = Q)
Inds_extremes<-which(Sum_found>Th)
THETA_found<-Estim_param_RF_homogeneous(Params_risk_Function = List_Params_RF,
                           Vect_l_function =MEV_FRECHET,
                           d=d,Q_thresh = Q)
THETA_found
print(THETA)
Gap<-abs(THETA_found$alpha-THETA)
print(paste0("The absolute difference between the estimator and the true value is ",Gap))
Exts_mev<-MEV_FRECHET[Inds_extremes,]
# if(sum(W_chosen)==1){
#   png(filename = paste0("functions_MV/Weighted_results_Lp=",q,".png"))
#   plot(MEV_FRECHET,log="xy",xlab="First dimension",ylab="Second dimension",
#        cex.lab=1.5)
#   abline(h=Th,v=Th,lty="dashed",col="blue")
# }else{
#   png(filename = paste0("functions_MV/results_Lp=",q,".png"))
#   plot(MEV_FRECHET,log="xy",xlab="First dimension",ylab="Second dimension",
#        cex.lab=1.5)
#   abline(h=Th,v=Th,lty="dashed",col="blue")
# }
## Change the threshold used.
M<-2000
#### Directly gives the lambda square coefficients for rgpapr != mev, gives the sigma...
# Simuls_gp<-Simulation_gParetoP(Params_risk_Function=List_Params_RF,
#                     Threshold= as.numeric(Th),d=ncol(MEV_FRECHET),
#                     M=M,theta_opt=THETA_found)
Simuls_gp<-mev::rgparp(n = M,shape = 1,thresh = as.numeric(Th),
            risk = "max",d = ncol(MEV_FRECHET),
            model = model_used,param = THETA,
            scale =rep(1,ncol(MEV_FRECHET)),
            loc=rep(1,ncol(MEV_FRECHET) ))
#png("test.png")
plot(Exts_mev,log="xy")
points(Simuls_gp,col="red")
data_forgg<-rbind.data.frame(Exts_mev,Simuls_gp)
data_forgg$origin<-c(rep("obs",nrow(Exts_mev)),
                     rep("simuls",nrow(Simuls_gp)))
require(ggside)
ggplot(data=data_forgg,aes(x=V1,y=V2,col=origin))+
  geom_point()+
  scale_x_continuous(transform="log10")+
  scale_y_continuous(transform="log10")+
  geom_xsidedensity(data=data_forgg,aes(fill=origin), 
                    alpha = 0.5)+
  geom_ysidedensity(data=data_forgg,aes(fill=origin), 
                    alpha = 0.5)

#dev.off()
# Unif_scale_sim<-apply(X = Simuls_gp,MARGIN=2,FUN = evd::pgpd,loc=0,
#       scale=1, shape=1)
### Scale of data observations
Unif_scale_exts<-exp(-(Exts_mev)^(-1))
### Enigma) what is the distrib of MEV ? 
#Unif_scale_sim<-exp(-Simuls_gp^(-1))
Unif_scale_sim<-apply(Simuls_gp,MARGIN=2,FUN = evd::pgpd,loc=0,
                     scale=1, shape=1)

plot(Unif_scale_exts,log="xy")
points(Unif_scale_sim,col="red",
       cex=0.5)
dev.off()
print("Kendall data vs simul")
print(cor(x = Exts_mev[,1],y=Exts_mev[,2],method = "kendall"))
print(cor(x=Simuls_gp[,1],y=Simuls_gp[,2],method="kendall"))

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
