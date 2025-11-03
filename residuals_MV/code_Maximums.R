# # # second_option -----------------------------------------------------------
# # # parametre clé !!
# # BEG<-10
# # END<-280
# # Seq_b<-seq.int(from = BEG,to = END,by = 5)
# # Result_AD_KS<-rbind.data.frame(sapply(Seq_b,AD_KS_B_choice,
# #                                       Vect_l_function = Vect_l_function,
# #                                       l_name=l_name))
# # Rtest<-as.data.frame(apply(X = Result_AD_KS,MARGIN = 1,FUN = unlist))
# # Rtest$variable<-l_name[as.numeric(Rtest$index_col)]
# # # unlist sur chaque colonne pr extraire la matrice
# # head(Rtest,10)
# # rownames(Rtest)<-c(1:nrow(Rtest))
# # Rtest$NPY_BLOCK<-as.numeric(Rtest$NPY_BLOCK)
# # Rtest$value<-as.numeric(Rtest$value)
# # Rtest$variable<-ifelse(Rtest$variable=="Surcote",
# #                        "S",Rtest$variable)
# # Rtest$index_col<-NULL
# # NB_YEARS<-37
# # Ratio_per_year<-round(nrow(Vect_l_function)/NB_YEARS)
# # ggplot(data=Rtest,aes(x=NPY_BLOCK,y=value,col=name_test))+
# #   facet_wrap(~variable)+
# #   geom_point()+
# #   geom_line()+
# #   xlab("Block length")+
# #   ylab("p value Frechet law")+
# #   geom_hline(yintercept=0.05)+
# #   labs(col="Legend")+
# #   geom_vline(xintercept=Ratio_per_year,
# #              col="orange",linetype="dashed")+
# #   annotate(x=Ratio_per_year,y=0.8,
# #            label="year block",geom="text")
# # #text pour ne pas avoir encadré, 
# # #label si on veut. 
 
# # mev::taildep(data = Vect_input,confint = "wald")
# # mev::taildep(data = Vect_l_function,confint = "wald")
# # 
# # # Compare with simulations. 
# # # # Max-stable model
# # # dat <- rmev(n = 5000, d = 2, param = 0.6, model = "log")
# # # goftest::ad.test(dat[,1],null = extRemes::"pevd",
# # #                  scale=1,loc=1,shape=1,type="GEV")
# # # goftest::ad.test(dat[,2],null = extRemes::"pevd",
# # #                  scale=1,loc=1,shape=1,type="GEV")
# # # taildep(dat, confint = 'wald')
# # 
# 
# 
# # OLD NON param -----------------------------------------------------------
# 
# Bound<-min(apply(X = Unifs_l,MARGIN = 1,FUN = max))
# if(ncol(Excedents_lprime)==2){
#   Mod_cop<-VineCopula::BiCopSelect(u1 = Unifs_l[,1],
#                                    u2=Unifs_l[,2])
#   Prop_l<-BiCopSim(N=M,family =Mod_cop$family,par = Mod_cop$par,
#                    par2 =Mod_cop$par2,
#                    check.pars = TRUE)
#   Max_<-apply(X = Prop_l,MARGIN = 1,FUN = max)
#   issue<-which(Max_<Bound)
#   Nb_issue<-length(issue)
#   Sim_Unifs_l<-Prop_l[-issue,]
#   while(Nb_issue>0){
#     Prop_l2<-BiCopSim(N=Nb_issue,family =Mod_cop$family,par = Mod_cop$par,
#                       par2 =Mod_cop$par2,
#                       check.pars = TRUE)
#     if(Nb_issue>1){
#       Max_<-apply(X = Prop_l2,MARGIN = 1,FUN = max)
#       issue<-which(Max_<Bound)
#       Nb_issue<-length(issue)
#       print(Nb_issue)
#       if(Nb_issue==0){
#         Sim_Unifs_l<-rbind(Sim_Unifs_l,Prop_l2)
#       } 
#       else{
#         Sim_Unifs_l<-rbind(Sim_Unifs_l,Prop_l2[-issue,])
#       }
#       
#     }
#     else{
#       Prop_l2<-c(Prop_l2)
#       Max_<-max(Prop_l2)
#       if(Max_>Bound){
#         Sim_Unifs_l<-rbind(Sim_Unifs_l,
#                            Prop_l2)
#         Nb_issue<-0
#       }
#     }
#   }
# }
# else{
#   ## construct structure matrix (mandatory for d>2)
#   Matrice_C<-function_Structure_Matrice(NB_dim = ncol(DF))
#   Mod_cop<-VineCopula::RVineCopSelect(Unifs_l)
#   Prop_l<-VineCopula::RVineSim(N=M,RVM = Mod_cop)
#   Max_<-apply(X = Prop_l,MARGIN = 1,FUN = max)
#   issue<-which(Max_<Bound)
#   Nb_issue<-length(issue)
#   Sim_Unifs_l<-Prop_l[-issue,]
#   while(Nb_issue>0){
#     Prop_l2<-VineCopula::RVineSim(N=Nb_issue,RVM = Mod_cop)
#     if(Nb_issue>1){
#       Max_<-apply(X = Prop_l2,MARGIN = 1,FUN = max)
#       issue<-which(Max_<Bound)
#       Nb_issue<-length(issue)
#       Sim_Unifs_l<-rbind(Sim_Unifs_l,Prop_l2[-issue,])
#     }
#     else{
#       Max_<-max(Prop_l2)
#       if(Max_>Bound){
#         Sim_Unifs_l<-rbind(Sim_Unifs_l,
#                            Prop_l2)
#         Nb_issue<-0
#       }
#     }
#   }
# }
# Sim_Unifs_l<-matrix(Sim_Unifs_l,ncol=ncol(Unifs_l),
#        nrow=length(Sim_Unifs_l)/ncol(Unifs_l),byrow = TRUE)
# plot(Sim_Unifs_l,main="Simulated uniform coordinates of the l couple")
# Sim_l<-Sim_Unifs_l
# for(j in c(1:ncol(Vect_l_function))){
#   Sim_l[,j]<-quantile(Excedents_lprime[,j],
#                       Sim_Unifs_l[,j])
# }
# }
