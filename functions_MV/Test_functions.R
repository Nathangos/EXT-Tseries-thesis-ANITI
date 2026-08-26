rm(list=ls())
set.seed(133)
require(mev)
require(reshape2)
require(ggplot2)
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

### Test Quantile regression
N_mev<-1000
N_sim<-500
MStable_sims<-mev::rmev(n = N_mev,d = 2,param = 0.7,model = "log")
Tau_seqce<-seq.int(from = 0.5,to = 0.95,
                   by = 0.1)
Exp_sims<-as.data.frame(-log(1-exp(-MStable_sims^(-1))))
colnames(Exp_sims)<-c("X_regressor","Y_target")
New_X<-rexp(N_sim)
BASISF_qreg<-"tp"
Qreg_lines<-lapply(Tau_seqce,
                   ALD_reg_XY,data_XY = Exp_sims,
                   basis_functionGAM =BASISF_qreg ,Sim_X=New_X
)
Combinations<-do.call(rbind,Qreg_lines)
Combinations$Tau<-as.numeric(Combinations$Tau)
list_ggplot<-list()
First_name<-"V1"
Scd_name<-"V2"
for(j in c(1:length(unique(Combinations$origin)))){
  cat<-unique(Combinations$origin)[j]
  ind_spe<-which(Combinations$origin==cat)
  Sub_comb<-Combinations[ind_spe,]
  Orig<-ggplot(data=Exp_sims,aes(x=X_regressor,
                                    y=Y_target))+
    geom_point()+
    xlab(First_name)+ylab(Scd_name)
  if(cat=="first"){
    Half<-Orig+
      geom_line(data=Sub_comb,aes(x=X_reg,y=Y_target,
                                  group=interaction(Tau,origin),
                                  col=Tau,linetype=origin))+
      guides(linetype="none")+
      scale_color_continuous(palette = "viridis")
    
    
  }else{
    Half<-Orig+
      coord_flip()+
      geom_line(data=Sub_comb,aes(x=X_reg,y=Y_target,
                                  group=interaction(Tau,origin),
                                  col=Tau,linetype=origin))+
      guides(linetype="none")+
      scale_color_continuous(palette = "viridis")
    
    
  }
  Full<-Half+
    geom_abline(slope = 1,intercept = 0,linetype="dashed")+
    theme(axis.title=element_text(size=20),
          legend.text=element_text(size=12),
          legend.title = element_text(size=13),
          axis.text = element_text(size=14))
  ggsave(filename = paste0(getwd(),"/functions_MV/test_MEV_qreg_",j,".png"),
         plot = Full,width=8,height=6)
}

#### Graphics to explain confidence bands for return level lines
levels_used<-c(0.7)
Curves_found<-sapply(levels_used,FUN = function(x){
  CURVE<-texmex::JointExceedanceCurve(Sample = Exp_sims,
                               ExceedanceProb = x)
  return(CURVE)
})
Curves_found_1sample<-do.call(cbind.data.frame,
  Curves_found)
colnames(Curves_found_1sample)<-colnames(Exp_sims)
Fction_new__fast_new_RL<-function(obs){
  Inds_t_s<-sample(x = c(1:nrow(obs)),size = nrow(obs),
                   replace = TRUE)
  New_data<-obs[Inds_t_s,]
  Curves_found_1sample<-sapply(levels_used,FUN = function(x){
    texmex::JointExceedanceCurve(Sample =  New_data,
                                 ExceedanceProb = x)
  })
  df<-as.data.frame(Curves_found_1sample)
  return(df)
}
N_test<-20
Ex_biv_rl_curves<-replicate(N_test,
            expr = Fction_new__fast_new_RL(Exp_sims))
GG_for_explan<-ggplot()+
  geom_jointExcCurve(x = Curves_found_1sample,
                     aes(x=X_regressor,y=Y_target,
                     col="estimator_data"))
palette_colors <- colors()
require(texmex)
for(z in c(1:N_test)){
  GG_for_explan<-GG_for_explan+
    geom_jointExcCurve(x = Ex_biv_rl_curves[[z]],
                aes(x=X_regressor,y=Y_target,
                 ), col=palette_colors[z],
                alpha=0.7)
}
theta_vect<-seq.int(from = 0,to = pi/2,length.out = 5)

for(theta in theta_vect){
  GG_for_explan<-GG_for_explan+
  geom_abline(slope =tan(theta),intercept = 0,col="black")
}
GG_for_explan<-GG_for_explan+
  labs(col="Legend")

ggsave(paste0(getwd(),"/functions_MV/showMethod_confRL.png"),
       width=8,height=6,plot = GG_for_explan)
