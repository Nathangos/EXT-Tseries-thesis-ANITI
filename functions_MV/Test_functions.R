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
Run_diagnostics_Gamma_G(LIMS_Y =LIMS_Y,
                        dims_elt_text = DIMS_plot,
                        Vectors_HTAIL = Simul_mev,
                        I = 1,J = 2,NAME_Vars = Name_vars,
                        Vect_k = Vect_k,q=0.8)
