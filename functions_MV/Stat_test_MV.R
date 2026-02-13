
# Test for exponent measure -----------------------------------------------
Stat_Test_Measure<-function(Boot_rep,Mod_BIV,Size_data,w,MAR1,
                            MAR2){
  vect_distances<-replicate(n=Boot_rep,Sample_stat_value(Size_data,Mod_BIV,w,
                                                         MAR1,
                                                         MAR2))
  return(vect_distances)
}

Sample_stat_value<-function(Size_data,Mod_BIV,w,MAR1,
                            MAR2){
  Simuls<-evd::rbvevd(n = Size_data,
                      asy=c(Mod_BIV$estimate[["asy1"]],
                            Mod_BIV$estimate[["asy2"]]),
                      dep = Mod_BIV$estimate[["dep"]],
                      model = Mod_BIV$model,
                      alpha=Mod_BIV$estimate[["alpha"]],
                      beta=Mod_BIV$estimate[["beta"]],
                      mar1 = MAR1,
                      mar2 = MAR2)
  ### Estimate model parameter Theta
  Estimated_model<-evd::fbvevd(x =  Simuls,
                               method="BFGS",model=Mod_BIV$model,
                               std.err = FALSE)
  New_Aparam<-evd::abvevd(w,
                          asy=c(Estimated_model$estimate[["asy1"]],
                                Estimated_model$estimate[["asy2"]]),
                          dep = Estimated_model$estimate[["dep"]],
                          model = Estimated_model$model,
                          alpha=Estimated_model$estimate[["alpha"]],
                          beta=Estimated_model$estimate[["beta"]])
  New_Anonparam<-evd::abvnonpar(w,data=Simuls,
                                method="pick")
  delta<-diff(w)[1]
  Sn_boot<-sum(delta*nrow(Simuls)*(New_Aparam-New_Anonparam)**(2))
  return(Sn_boot)
}
