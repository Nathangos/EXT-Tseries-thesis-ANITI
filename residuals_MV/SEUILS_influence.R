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
source("fonctions/fonctions_Pareto.R")
#source("Modele_ss_tend/Travail_angles.R")
source("fonctions/fonctions.r")
source("fonctions/fonctions_r_pareto_EVOLTAU.r")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

Influence_Threshs<-function(name_variable,type_entree,coeurs,lien_racine,
                            Q1_th,Q2_th,NT_th,N_v){
  lien_donnees<-paste0(lien_racine,name_variable,"_residus.csv")
  if(type_entree=="clust"){
    All<-read.csv(file=lien_donnees)[,2:39]
    colnames(All)<-c(1:38)
    Vrais_indices<-All[,38]
    Donnes<-All[,c(1:37)]
  }
  else{
    Donnes<-read.csv(file=lien_donnees)[,2:38]
  }
  colnames(Donnes)<-c(1:37)
  rownames(Donnes)<-c(1:nrow(Donnes))
  Result_all_times<-list()
  for(t in c(1:ncol(Donnes))){
   Result_all_times[[t]]<-tryCatch(Outils_Threshr_choix(Q1 = Q1_th,
                                                        Q2 = Q2_th,
                                 NT_ths = NT_th,
                                 variable = Donnes[,t],
                                 N_v=N_v)
            ,warning=function(w) return(
              list(Outils_Threshr_choix(Q1 = Q1_th,Q2 = Q2_th,
                                        NT_ths = NT_th,
                                        variable = Donnes[,t],N_v=N_v),w)))
   
  }
  return(Result_all_times)
}


# Parametres --------------------------------------------------------------
l_root<-"residus_clust/"
Q1<-0.70
Q2<-0.95
Number_candidates<-20
TE<-"residus"
Quantiles_candidates<-seq.int(Q1,Q2,length.out = Number_candidates)
NOMV<- "Hs"

# Lancement ---------------------------------------------------------------
Result_Analysis_BiGPD<-Influence_Threshs(name_variable =NOMV,coeurs = coeurs,
                  lien_racine = l_root,Q1_th = Q1,
                  type_entree = TE,
                  Q2_th = Q2,NT_th = Number_candidates,N_v=1)
vect_w<-rep(NA,37)
Mean_w<-rep(NA,37)
par(mfrow=c(2,2))
pbs<-rep(NA,37)
for(t in c(1:37)){
  r_t<-Result_Analysis_BiGPD[[t]]
  if(length(r_t)==2){
    x<-plot(r_t[[1]])$x
    y<-plot(r_t[[1]])$y
    z<-y/sum(y)
    vect_w[t]<-which.max(y)
    Mean_w[t]<-sum(x*z)/100
    pbs[t]<-"pb_conv"
  }
  else{
    x<-plot(r_t)$x
    y<-plot(r_t)$y
    z<-y/sum(y)
    Mean_w[t]<-sum(x*z)/100
    vect_w[t]<-which.max(plot(r_t)$y)
  }
}
Q_chosen<-Quantiles_candidates[vect_w]
Q_chosen
par(mfrow=c(1,1))

# Export of the thresholds ------------------------------------------------
write.csv(x=Mean_w,file = paste0("residus_MV/seuil_",NOMV,".csv"))

# End of the code ---------------------------------------------------------
# -- ----------------------------------------------------------------------

stopCluster(coeurs)
