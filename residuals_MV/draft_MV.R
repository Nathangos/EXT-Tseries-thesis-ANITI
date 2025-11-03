# EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")
# db<-as.data.frame(Excedents_lprime[,Z]-1)
# colnames(db)<-c("x")
# con <- gamlss.control(n.cyc = 100,
#                       mu.step = 0.1, 
#                       sigma.step = 0.1, nu.step = 0.1,tau.step = 0.1,autostep=TRUE)
# con.i=glim.control(glm.trace = FALSE)

# # Identity to enable mu/gamma<0
# Fitting_riskF_z<- gamlss(x~1, 
#                          data=db,
#                          family = EGPD1Family(mu.link = "identity"),
#                          control = con,mu.start=INIT[1],sigma.start=INIT[2],
#                          nu.start=0.5,
#                          i.control=con.i,
#                          method=CG())
# muFit <- fitted(Fitting_riskF_z,"mu")[1]
# sigmaFit <- predict(Fitting_riskF_z,what="sigma", 
#                     type="response")[[1]]
# nuFit <- predict(Fitting_riskF_z,what="nu", 
#                  type="response")[[1]]


# for(Z in c(1:2)){
#   R_X<-sapply(Excedents_lprime[,Z],CDF_1d_Mbiv,z=Z,
#               Params_Biv=Params_Biv,Th=Seuil_lprime)
#   Index<-which(is.na(R_X)==TRUE)
#   Subset<-R_X
#   if(length(Index)>0){
#     Subset<-R_X[-Index]
#   }
#   plot(density(Subset))
#   print(goftest::ad.test(Subset,null = "punif")$p.value)
# }
# 
# 
# if(Model!="hr"){
#   Sim_l<-mev::rgparp(n = M,shape = 1,
#                      thresh = Seuil_lprime,risk = "max",
#                      model = Model,d = ncol(Excedents_lprime),
#                      param = Alpha_opt,
#                      scale = scale_frechet_d,
#                      loc = scale_frechet_d)
# }else{
#   Sigma<-matrix(0,nrow = ncol(Excedents_lprime),
#                 ncol=ncol(Excedents_lprime))
#   Sigma[1,2]<-(1/Alpha_opt)^(2)
#   Sigma[2,1]<-Sigma[1,2]
#   Sim_l<-mev::rgparp(n = M,shape = 1,
#                      thresh = Seuil_lprime,risk = "max",
#                      model = Model,d = ncol(Excedents_lprime),
#                      param = Alpha_opt,
#                      scale = scale_frechet_d,
#                      loc = scale_frechet_d,
#                      sigma=Sigma)


# P_BIV<-list("model"=Mod_biv$model,
#             "cste"=scale_frechet_d,
#             "dep"=ALPHA_Mstable)

# if(typemod%in%c("alog","aneglog")){
#   asy1<-Mod_biv$estimate[["asy1"]]
#   if(SYMbool==TRUE){
#     asy2<-asy1
#   }
#   if(SYMbool==FALSE){
#     asy2<-Mod_biv$estimate[["asy2"]]
#   }
#   ASY<-c(asy1,asy2)
# }
# if(typemod=="log"){
#   ASY<-rep(1,2)
# }
# P_BIV<-NULL
# if(TML=="param"){
#   z<-1
#   if("dep"%in%names(Mod_biv$estimate)){
#     Dep<-Mod_biv$estimate[["dep"]]
#   }
#   else{
#     Dep<-NA
#   }
#   Aparam<-evd::abvevd(w,dep = Dep,
#                       asy = ASY,
#                       beta = Mod_biv$estimate[["beta"]],
#                       alpha = Mod_biv$estimate[["alpha"]],
#                       model = Mod_biv$model)
#   Sn_Distance<-sum(nrow(Result_std)*(Aparam-Anonparam)**(2)*delta)
#   Sn_Distance
# 
#   MAR1<-c(Mod_biv$estimate[["loc1"]],
#           Mod_biv$estimate[["scale1"]],
#           Mod_biv$estimate[["shape1"]])
#   MAR2<-c(Mod_biv$estimate[["loc2"]],
#           Mod_biv$estimate[["scale2"]],
#           Mod_biv$estimate[["shape2"]])
#   source("residus_MV/Stat_test_MV.R")
#   Boot_rep<-200
#   Distrib_Sn_pb<-Stat_Test_Measure(Boot_rep = Boot_rep,Mod_BIV = Mod_biv,
#                                    Size_data = nrow(Result_std),
#                                    SYMbool = SYMbool,
#                                    w = w,MAR1=MAR1,MAR2=MAR2)
#   summary(Distrib_Sn_pb)
#   plot(density(Distrib_Sn_pb))
#   abline(v=Sn_Distance,col="red")

# }
# W_chosen<-c(0.7,0.3)
# Risk_f<-function(x){return(sum_penalized(x = x,vect_w = W_chosen))}
# Name_riskF<-"sum_penalized"
# wFun<-function(x,u){
#   weightFun_sum_penalized(x = x,u = u,
#                           vect_w = W_chosen)
# }
# dwFun<-function(x,u){
#   DweightFun_sum_penalized(x = x,u = u,
#                            vect_w = W_chosen)
# }

Risk_f<-function(x){return(sum(x))}
Name_riskF<-"sum"
wFun<-weightFun_sum
dwFun<-DweightFun_sum

scale_frechet_d<-c(1,1)
List_Params_RF<-list("parametric"=TRUE,"name_RF"=Name_riskF,
                     "function"=Risk_f,"name_model"="log",
                     "Scale_Frechet"=scale_frechet_d,
                     "weightFunction"=wFun,
                     "dWeightFunction"=dwFun,
                     "init_opt_param"=0.5)
Theta_chosen<-0.7
source("functions_MV/function_estim_simul_gPareto.R")
simul_g<-Simulation_gParetoP(Params_risk_Function = List_Params_RF,
                    Threshold = 20,d = 2,
                    M = 10000,theta_opt = Theta_chosen)
par(mfrow=c(1,2))
POT::chimeas(simul_g,ask=FALSE)
par(mfrow=c(1,1))
simul_g<-as.data.frame(simul_g)
plot(simul_g,log="xy")
df<-as.data.frame(simul_g[,1])
colnames(df)<-c("x")
Model_gp<-mev::gp.fit(xdat = simul_g,threshold = 0)

# POT property, fractal  --------------------------------------------------
RF_data<-apply(X = simul_g,MARGIN = 1,FUN =Risk_f)
u_lg<-as.numeric(quantile(RF_data,probs =0.95))
Inds_exts<-which(RF_data>u_lg)
POT_taken<-simul_g[Inds_exts,]
plot(POT_taken,log="xy")
Theta_opt<-Estim_param_RF_homogeneous(Params_risk_Function = List_Params_RF,
                           Seuil_lprime = u_lg,Vect_l_function = simul_g)
Theta_opt
simul_estim<-Simulation_gParetoP(Params_risk_Function = List_Params_RF,
                    Threshold = u_lg,d = 2,
                    M = 100,theta_opt = Theta_opt)
points(simul_estim,col="red"
       ,cex=0.5)

# problem with the gradient scoring rule ----------------------------------
# Use result from resnick -------------------------------------------------
Ntest<-1000
GAMMA_0<-1
B<-rbinom(n =Ntest,size = 1,prob = 0.5)
Theta_1<-cbind(1+rexp(n = Ntest,rate = 1),rep(1,Ntest))
Theta_2<-cbind(rep(1,Ntest),1+rexp(n = Ntest,rate = 1))
Theta_null<-B*Theta_1+(1-B)*Theta_2
R_Pareto<-(runif(n = Ntest))^(-GAMMA_0)
V<-R_Pareto*Theta_null

# Multivariate regularly varying ------------------------------------------
plot(V,log="xy")
RF_data<-apply(X = V,MARGIN = 1,FUN = Risk_f)
ug<-quantile(RF_data,probs = 0.75)
Inds_exts<-which(RF_data>ug)
plot(V[Inds_exts,],log="xy")
par(mfrow=c(1,2))
POT::chimeas(V,ask=FALSE)
par(mfrow=c(1,1))

goftest::ad.test(x = V[,2],null = extRemes::"pevd",
                 threshold=1,shape=1,type="GP",scale=1)$p.value
## Not pareto
Fit_resnick<-mev::gp.fit(xdat = (V[,1]),threshold = 0)
Theta_found<-Fit_resnick$est
mev::fit.extgp(data = (V[,1]),model=1,method="pwm",
               init=c(0.5,Theta_found))
Graphics_estimators_gamma(series = V[,1],
                          vect_k = c(10:300),
                          Title_graphic = "test Resnick shape")
Theta_opt_Resnick<-Estim_param_RF_homogeneous(Params_risk_Function = List_Params_RF,
                           Seuil_lprime = ug,Vect_l_function = V)
Theta_opt_Resnick
