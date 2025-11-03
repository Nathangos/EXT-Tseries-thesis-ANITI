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
### Source des fonctions. ####
#Functions from univariate case. 
source("../fonctions/fonctions_Pareto.R")
source("../fonctions/fonctions.r")
source("../fonctions/fonctions_r_pareto_EVOLTAU.r")

### Generation functions
source("functions_MV/function_generation_MV.R")
### Source des fonctions. ####
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)


G_used<-"Weighted_NQ"
q<-20
G_used2<-G_used
if(G_used=="Weighted_NQ"){
  G_used2<-paste0(G_used2,"_q=",q)
}

Link_resid<-paste0("residuals_MV/")
lNAME<-c("U","Surcote")
distance_sq<-function(x){
  return(sum(x^2)**(1/2))
}
link_export_<-paste0("Generated_extremes/",G_used2,"/")
try({
  dir.create(file.path(link_export_))
})
cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue",
         "model"="red","Storm"="black", "exceedances"="red")
Output<-Generator_TSERIES_MV(type_donnees = "HIVER",coeurs = coeurs, 
            list_variable = lNAME,
            link_import_residuals = Link_resid,
            link_export_generations=link_export_,
            link_VAR_model = "whitening_MV/VAR_model.csv",
            K_chosen = 20,f_distance = distance_sq,
            cols_gg = cols_,Name_riskF=G_used2)

for(Nvariable in names(Output[[1]])){
  par(mfrow=c(1,2))
  EPSI<-Output$Epsi_sim
  X0<-Output$X0
  Target<-Output$Target
  Target_nv<-Target[[Nvariable]]
  Epsi_nv<-EPSI[[Nvariable]]
  max_y<-max(max(apply(X =Target_nv,MARGIN = 2,FUN = max)),
             max(apply(X =Epsi_nv,MARGIN = 2,FUN = max))
         )
  min_y<-min(min(apply(X =Target_nv,MARGIN = 2,FUN = min)),
             min(apply(X =Epsi_nv,MARGIN = 2,FUN = min)))
  matplot(t(Epsi_nv),type="l",
          ylim=c(min_y,max_y))
  matplot(t(Target_nv),type="l",
          ylim=c(min_y,max_y))
  mtext("Epsi sim and X0",line=-50)
  mtext("X sim vs ext X",line=-50)
  par(mfrow=c(1,1))
}

Data_U<-Output$Orig$Surcote
summary(Data_U)
Inds_exts<-Output$Inds_extremes$indexes
summary(Data_U[Inds_exts,])
### Comparing sim with obs

