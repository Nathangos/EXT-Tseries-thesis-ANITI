rm(list=ls())
source("residus_MV/F_simulations_cond.R")
set.seed(133)

Th<-40
DEP<-c(0.8,0.7)
MODEL<-"log"
Scale_frech<-c(0.98,1.12)
P_BIV<-list("model"=MODEL,
            "cste"=Scale_frech,
            "dep"=DEP)
Simuls_<-mev::rgparp(n = 1000,
            shape = 1,
            thresh = Th,
            risk = "max",
            d=2,loc = Scale_frech,
            scale = Scale_frech,
            param=DEP[1],
            model = MODEL)
plot(Simuls_,log="xy")
abline(h=Th,col="red",v=Th)
Z<-1
R_X<-sapply(Simuls_[,Z],CDF_1d_Mbiv,z=Z,
            Params_Biv=P_BIV,Th=Th)
summary(R_X)
Index<-which(is.na(R_X)==TRUE)
Simuls_[Index,Z]
Subset<-R_X
if(length(Index)>0){
  Subset<-R_X[-Index]
}
plot(density(Subset))
goftest::ad.test(Subset,null = "punif")

