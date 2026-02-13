# Other methods with transf and then heavy tail with Naveau code.

# PARETO<-do.call(cbind.data.frame,Resultat_P$obs)
# colnames(PARETO)<-c(1:ncol(PARETO))
# CHOSEN_XI<-(1/2)
# TRANSFO<-log(PARETO**CHOSEN_XI)
# SIGMA<-1
# INIT_GPD<-c(SIGMA,CHOSEN_XI)
# for(t in c(1:length(Resultat_P$gamma))){
#   VAR_t<-TRANSFO[,t]
#   Model_EXTGP<-mev::fit.extgp(data =VAR_t,
#                               model = 1,
#                               method = "mle",
#                               init =c(0.5,INIT_GPD),
#                               plots = FALSE)
#   Theta_extgp_t<-Model_EXTGP$fit$mle
#   
#   # Take estimation ---------------------------------------------------------
#   UNIF_t<-sapply(X =VAR_t,
#                  FUN=function(x){
#                    mev::pextgp(q =x ,kappa =Theta_extgp_t[["kappa"]],
#                                xi =Theta_extgp_t[["xi"]] ,
#                                sigma = Theta_extgp_t[["sigma"]],
#                                type = 1)})
#   MATRIX_TRANSFO[,t]<-UNIF_t
#   
# }