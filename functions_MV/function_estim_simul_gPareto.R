
Psi_gamma<-function(Mat_gamma,vect_ij){
  i<-vect_ij[1]
  j<-vect_ij[2]
  ref_j<-Mat_gamma[j,1]
  ref_i<-Mat_gamma[i,1]
  num_ij<-Mat_gamma[i,j]
  return((1/4*(ref_i+ref_j-num_ij)))
}

GAMMA_to_SigmaTilde<-function(Mat_GAMMA,k_const){
  MATRIX_nul<-matrix(0,nrow = nrow(Mat_GAMMA),
                     ncol=ncol(Mat_GAMMA))
  Coords<-t(utils::combn(x = c(1:ncol(Mat_GAMMA)),m = 2))
  Diag_coords<-t(sapply(c(1:ncol(Mat_GAMMA)),
                        FUN=function(x){return(rep(x,2))}))
  Coords<-rbind(Coords,Diag_coords)
  Sigma_values<-apply(X = Coords,MARGIN = 1,
                      Wrapper_function,d=ncol(MATRIX_nul),
                      Gam=Mat_GAMMA,k_taken=k_const)
  MATRIX_nul[Coords]<-Sigma_values
  MATRIX_nul<-MATRIX_nul+t(MATRIX_nul)
  diag(MATRIX_nul)<-diag(MATRIX_nul)/2
  Sigma_tilde<-MATRIX_nul
  return(Sigma_tilde)
}
# Simulation of pareto process when g=max ---------------------------------
Max_Likely<-function(alpha_value,scale_frechet_d,obs,
                     Model, 
                     threshold){
  Likely_data<-likmgp(
    obs,
    threshold,
    scale_frechet_d,
    scale_frechet_d,
    c(1,1),
    list("alpha"=alpha_value),
    model = Model,
    likt = "mgp",
    lambdau = 1
  )
  return(-Likely_data)
}

Simul_RF_max_NP<-function(Params_risk_Function,M,Seuil_lprime,Vect_l_function){
  Lg<-apply(X =Vect_l_function,MARGIN = 1,FUN = max)
  Inds_exts_lg<-which(Lg>Seuil_lprime)
  Excedents_lprime<-Vect_l_function[Inds_exts_lg,]
  Excedents_Exp_scale<-Vect_l_function
  Null_hyp_result<-matrix(NA,nrow = 3,
                          ncol = ncol(Excedents_lprime))
  MOD_EGPD_RISKF<-list()
  for(Z in c(1:ncol(Excedents_lprime))){
    INIT<-c(Seuil_lprime,1)
    Fitting_riskF_z<-mev::fit.extgp(data = Vect_l_function[,Z]-1,
                                    model = 1,
                                    method = "mle",
                                    init = c(0.5,INIT))
    Theta_extgpd<-Fitting_riskF_z$fit$mle
    sigmaFit<-Theta_extgpd[["sigma"]]
    nuFit<-Theta_extgpd[["kappa"]]
    muFit<-Theta_extgpd[["xi"]]
    Seuil_lprime_unif<-pEGPDModel1(q = Seuil_lprime-1,mu =muFit,
                                   sigma = sigmaFit,
                                   nu = nuFit)
    Unif_RiskF_z<-pEGPDModel1(q = Vect_l_function[,Z]-1,
                              mu =muFit,sigma = sigmaFit,
                              nu = nuFit)
    plot(density(Unif_RiskF_z),main=paste0("Dens for the ",Z,"th function"))
    MOD_EGPD_RISKF[[Z]]<-list("mu"=muFit,
                              "sigma"=sigmaFit,"nu"=nuFit)
    Excedents_Exp_scale[,Z]<--log(1-Unif_RiskF_z)
    Seuil_lprime_exp<--log(1-Seuil_lprime_unif)
    KS_<-ks.test(x = Excedents_Exp_scale[,Z],
                 "pexp",rate=1)$p.value
    CVM_<-goftest::cvm.test(Excedents_Exp_scale[,Z],"pexp",
                            rate=1)$p.value
    AD_<-goftest::ad.test(Excedents_Exp_scale[,Z],null="pexp",
                          rate=1)$p.value
    Null_hyp_result[,Z]<-c(KS_,CVM_,AD_)
  }
  # Simulations with LEGRAND method -------------------------------------------------------
  Sim_l<-Non_param_LEGRAND_RiskF(Data_scale_exp =Excedents_Exp_scale[Indices_exts,],
                                   nb_simul = M,MOD_EGPD_RiskF = MOD_EGPD_RISKF,
                                   Threshold_EXP=Seuil_lprime_exp)
    return(Sim_l)
}

#' Non_param_LEGRAND_RiskF
#'
#' @param Data_scale_exp: dataframe. Vector of risk functions for extreme individuals.
#' @param nb_simul: int. Number of desired simulations.  
#' @param MOD_EGPD_RiskF: list. EGPD parameters for each variable.
#' @param Threshold_EXP: float. Threshold used to identify exponential scale.
#'
#' @return Simulations of Risk functions. 
#' @export
#'
#' @examples
Non_param_LEGRAND_RiskF<-function(Data_scale_exp,nb_simul,
                                  MOD_EGPD_RiskF,Threshold_EXP){
  # From algorithm 1 in LEGRAND ---------------------------------------------
  # Simul Z exp -------------------------------------------------------------
  Zsimul<-rexp(n = nb_simul,rate = 1)+as.numeric(Threshold_EXP)
  D<-ncol(Data_scale_exp)
  
  # Simul T  ----------------------------------------------------------------
  Simul_Z<-matrix(rep(Zsimul,D),ncol = D,
                  nrow=nb_simul,
                  byrow = FALSE)
  
  Delta_i<-as.numeric(t(diff(t(Data_scale_exp))))
  
  # Sample Delta ------------------------------------------------------------
  Delta_tilde<-sample(Delta_i,size = nb_simul,
                      replace = TRUE)
  IND_plus<-which(Delta_tilde>0)
  Delta_mat<-matrix(0,ncol = D,
                    nrow=nb_simul)
  Delta_mat[-IND_plus,1]<-Delta_tilde[-IND_plus]
  Delta_mat[IND_plus,2]<-(-1)*Delta_tilde[IND_plus]
  
  # Simul complement --------------------------------------------------------
  Simul_RiskF<-Simul_Z+Delta_mat
  print(summary(Simul_RiskF))
  
  # Exp to orig -------------------------------------------------------------
  Simul_Z_Unif<-1-exp(-Simul_RiskF)
  print(summary(Simul_Z_Unif))
  Simul_RiskF<-Simul_Z_Unif
  for(Z in c(1:ncol(Simul_RiskF))){
    Orig_z<-Simul_Z_Unif[,Z]
    Theta_EXTGP<-MOD_EGPD_RiskF[[Z]]
    Orig_to_EGPD<-qEGPDModel1(p = Orig_z,mu = Theta_EXTGP[["mu"]],
                              sigma = Theta_EXTGP[["sigma"]],
                              nu = Theta_EXTGP[["nu"]])
    Simul_RiskF[,Z]<-Orig_to_EGPD
  }
  return(Simul_RiskF)
}

#' Title
#'
#' @param Params_risk_Function 
#' @param M 
#' @param Seuil_lprime 
#' @param Vect_l_function 
#'
#' @return Estimated parameter for the g Pareto law
#' @export
#'
#' @examples
Estim_param_RF_homogeneous<-function(Params_risk_Function,Seuil_lprime,Vect_l_function,
                                     d){
  Lg<-apply(X =Vect_l_function,MARGIN = 1,FUN = Params_risk_Function[["function"]])
  Inds_exts_lg<-which(Lg>Seuil_lprime)
  #as.list(as.data.frame(t(Exts_mev)))
  Excedents_lprime<-Vect_l_function[Inds_exts_lg,]
  if(Params_risk_Function[["name_RF"]]=="max"){
    Model<-Params_risk_Function[["name_model"]]
    scale_frechet_d<-Params_risk_Function[["Scale_Frechet"]]
    #Excedents_lprime<-sapply(X = Excedents_lprime,FUN = unlist)
    Optimization_<-optim(par =Params_risk_Function[["init_opt_param"]],
                         fn = Max_Likely,lower=0,upper = 1,
                         scale_frechet_d=scale_frechet_d,method = "Brent",
                        obs=Excedents_lprime,threshold=Seuil_lprime,
                         Model=Model)
    Alpha_Likely<-Optimization_$par
    return(Alpha_Likely)
    
  }else{
    # Set up, use of g and definition of the weight function  --------------------------
    Model<-Params_risk_Function[["name_model"]]
    #scale_frechet_d<-Params_risk_Function[["Scale_Frechet"]]
    Excedents_lprime<-as.list(as.data.frame(t(Excedents_lprime)))
    weightFun<-Params_risk_Function[["weightFunction"]]
    DweightFun<-Params_risk_Function[["dWeightFunction"]]
    
    # Optimisation of the objective function --------------------------------------------
    if(Model=="log"){
      #method = "L-BFGS-B",lower=0,upper=1
      Optimization_Score<-optim(par = Params_risk_Function[["init_opt_param"]],
                                fn = Objective_function,weightFun = weightFun ,
                                lower=0,upper=1,
                                method = "Brent",dWeightFun = DweightFun,Extreme_inds  = Excedents_lprime,
                                u=Seuil_lprime,
                                model_used=Model,
                                d=d)
      
    }
    else{
      THETA_init<-Params_risk_Function[["init_opt_param"]]
      v <- (d * (d - 1)) / 2
      # l<-1
      # lower_bounds <- rep(0, l + v)
      # larger than 0 as we have the squared lambda !
      # not sufficient to impose what we want...
      Optimization_Score<-optim(par =THETA_init ,
                                method = "L-BFGS-B" ,fn = Objective_function,
                                weightFun = weightFun,dWeightFun = DweightFun,
                                lower=rep(0,v),Extreme_inds  = Excedents_lprime,
                                u=Seuil_lprime,model_used=Model,
                                control = list(maxit = 1000),
                                d=d)
    }
    ALPHA_opt_scoring<-Optimization_Score$par
    return(ALPHA_opt_scoring) 
  }
}

#' Objective_function
#'
#' @param theta 
#' @param weightFun 
#' @param dWeightFun 
#' @param Extreme_inds 
#' @param u 
#' @param model_used 
#'
#' @return
#' @export
#'
#' @examples
Objective_function<-function(theta,weightFun,dWeightFun,
                             Extreme_inds,u,model_used,d){
  if(typeof(Extreme_inds)=="list"){
    n<-length(Extreme_inds)
  }
  if(typeof(Extreme_inds)=="double"){
    n<-nrow(Extreme_inds)
  }
  if(model_used=="hr"){
    # Take the gamma parameters
    # and fill the corresponding matrix Mat_Gamma
    # d parameters
    GAMMAnon_diag<-theta
    print(theta)
    Mat_GAMMA<-matrix(0,nrow = d,ncol=d)
    Coords<-t(utils::combn(x = c(1:d),m = 2))
    Mat_GAMMA[Coords]<-GAMMAnon_diag
    Mat_GAMMA<-t(Mat_GAMMA)+Mat_GAMMA
    diag(Mat_GAMMA)<-rep(0,ncol(Mat_GAMMA))
    # Sigma used in Hreiss Sigma.
    Sigma_<-Mat_GAMMA
    Diag_coords<-t(sapply(c(1:d),
                          FUN=function(x){return(rep(x,2))}))
    WHOLE_coord<-rbind(Coords, Diag_coords)
    # Sigma_<-GAMMA_to_SigmaTilde(Mat_GAMMA = Mat_GAMMA)
    Sigma_[WHOLE_coord]<-apply(WHOLE_coord,FUN = Psi_gamma,
                                   Mat_gamma=Mat_GAMMA,MARGIN = 1)
    Sigma_<-t(Sigma_)+Sigma_
    diag(Sigma_)<-diag(Sigma_)/2
    sigmaInv <- MASS::ginv(Sigma_)
    sigma <- diag(Sigma_)
    q <- rowSums(sigmaInv)
    A <- sigmaInv - q %*% t(q)/sum(q)
    zeroDiagA <- A
    diag(zeroDiagA) <- 0
    mtp <- 2 * q/(sum(q)) + 2 + sigmaInv %*% sigma - (q %*% t(q) %*% 
                                                        sigma)/(sum(q))
    theta<-list("A"=A,"mtp"=mtp)
    
  }
  scores <- lapply(1:n, computeScores_Gpareto,
                   model_used=model_used,theta=theta,Extreme_inds =Extreme_inds,
                   dWeightFun=dWeightFun,weightFun=weightFun,u=u)
  Result<-sum(unlist(scores))/n
  return(Result)
}

#' computeScores_Gpareto
#'
#' @param theta 
#' @param i 
#' @param model_used 
#' @param Extreme_inds
#' @param u
#'
#' @return
#' @export
#'
#' @examples
computeScores_Gpareto <- function(theta,i,model_used,Extreme_inds,
                                  weightFun,dWeightFun,d,u){
  if(typeof(Extreme_inds)=="list"){
    PTO_ind<-Extreme_inds[[i]]
  }
  if(typeof(Extreme_inds)=="double"){
    PTO_ind<-Extreme_inds[i,]
  }
  weights<-weightFun(x = PTO_ind,u = u)
  dWeights<-dWeightFun(x = PTO_ind,u = u)
  if(model_used=="log"){
    ALPHA<-theta
    DENOMTOR<-sum(PTO_ind^(-1/ALPHA))
    gradient<-(-1-(1/ALPHA))*PTO_ind^(-1)+
      (-(ALPHA-2)/ALPHA)*((PTO_ind^(-(1/ALPHA)-1))/DENOMTOR)
    NUMTOR1<-DENOMTOR*(-(1/ALPHA)-1)*PTO_ind^(-(1/ALPHA)-2)
    NUMTOR2<-(-1/ALPHA)*PTO_ind^((-(1/ALPHA)-1)*2)
    NUMTOR<-NUMTOR1-NUMTOR2
    diagHessian <- (1+(1/ALPHA))*PTO_ind^(-2)-
      ((ALPHA-2)/ALPHA)*(NUMTOR/DENOMTOR^(2))
    
    return(sum(2 * (weights * dWeights) * gradient + weights^2 * 
                 diagHessian + 1/2 * weights^2 * gradient^2))
  }
  if(model_used=="hr"){
    A<-theta[["A"]]
    mtp<-theta[["mtp"]]
    gradient <- -1/2 * ((A + t(A)) %*% log(PTO_ind)) * (1/PTO_ind) - 
      1/2 * (1/PTO_ind) * mtp
    diagHessian <- -1/2 * diag(A + t(A)) * (1/PTO_ind^2) + 
      1/2 * ((A + t(A)) %*% log(PTO_ind)) * (1/PTO_ind)^2 + 
      1/2 * (1/PTO_ind)^2 * mtp
  }
  return(sum(2 * (weights * dWeights) * gradient + weights^2 * 
               diagHessian + 1/2 * weights^2 * gradient^2))
}

#' Simulation_gParetoP
#'
#' @param Params_risk_Function: list. Parameters of the simulation.  
#' @param Threshold: float. Threshold.
#' @param d: int. Number of dimensions. 
#' @param theta_opt: float. Estimated estimator of the copula in
#' the chosen model. 
#'
#' @return Simulations of g pareto process for homogeneous functions.
#' @export
#'
#' @examples
Simulation_gParetoP<-function(Params_risk_Function,Threshold,d,M,
                              theta_opt){
  Nameg<-Params_risk_Function[["name_RF"]]
  name_model<-Params_risk_Function[["name_model"]]
  if(Nameg%in%c("sum","max","min","l2")){
    if(name_model=="log"){
      Simul_g<-mev::rgparp(n = M,shape = 1,risk =Nameg ,
                           d = d,
                           thresh=Threshold,
                           model = name_model,
                           param =theta_opt,
                           loc=rep(1,d),
                           scale=rep(1,d))
    }
    if(name_model=="hr"){
      Simul_g<-mev::rgparp(n = M,shape = 1,risk =Nameg ,
                           d = d,
                           thresh=Threshold,
                           model = name_model,
                           sigma = theta_opt,
                           loc=rep(1,d),
                           scale=rep(1,d))
    }
    
    return(Simul_g)
  }
  else{
    Simul_g<-replicate(M,Rejection_sampling_gPto(Params_risk_Function=Params_risk_Function,
                                        Threshold=Threshold,
                                        d=d,
                                        theta_opt=theta_opt))
    Simul_g<-matrix(Simul_g,ncol=d,byrow = TRUE)
  }
  return(Simul_g)
}

#' Rejection_sampling_gPto
#'
#' @param Params_risk_Function: list. Parameters of the simulation.  
#' @param Threshold: float. Threshold.
#' @param d: int. Number of dimensions. 
#' @param theta_opt: float. Estimated estimator of the copula in
#' the chosen model. 
#'
#' @return Simulations of g pareto process for functions not handled by the package
#' @export
#'
#' @examples
Rejection_sampling_gPto<-function(Params_risk_Function,Threshold,d,
                                  theta_opt){
  
  kept<-FALSE
  g_chosen<-Params_risk_Function[["function"]]
  Model<-Params_risk_Function[["name_model"]]
  while(kept==FALSE){
    unif_sample<-runif(1)
    Copy<-Params_risk_Function
    # U generator
    if(Model!="hr"){
      candidate<-mev::rmev(n =1,d = d,param = 1/theta_opt,
                       model = Model)
    }
    else{
      candidate<-mev::rmev(n=1,d = d,sigma = theta_opt,
                       model = Model)
    }
    # candidate<-Simulation_gParetoP(Params_risk_Function =Copy ,
    #                                Threshold = 1,
    #                                d = d,M = 1,
    #                                theta_opt = theta_opt)
    
    Cond_factor<-g_chosen(candidate)
    if(unif_sample<=Cond_factor){
      Angular_sim<-(candidate/Cond_factor)
      sim_pto<-(1-runif(1))^(-1)*Threshold
      Full_gpto<-(Angular_sim*sim_pto)
      return(Full_gpto)
    }
  }
}

# Weight functions for several risk functions g ---------------------------


# Sum ---------------------------------------------------------------------
weightFun_sum<-function(x,u){
  return(x * (1 - exp(1-sum(x/u))))
}
DweightFun_sum<-function(x,u){
  return((1 - exp(1-sum(x/u)))+ (x/u)*exp(1-sum(x/u)))
}
sum_penalized<-function(x,vect_w){
  return(sum(vect_w*x))
}
Shape_boundary<-function(x_j,j,ul,RiskF){
  xj_bound<-uniroot(f = F_foruniroot_boundary_thresh,
                    interval = c(-100,100),
                    x_j=x_j,j=j,RiskF=RiskF,ul=ul)$root
  return(xj_bound)
}
F_foruniroot_boundary_thresh<-function(xj,x_j,j,ul,
                                              RiskF){
  vect_xinput<-rep(NA,length(x_j)+1)
  vect_xinput[-j]<-x_j
  vect_xinput[j]<-xj
  opt_cost<-RiskF(vect_xinput)-ul
  return(opt_cost)
}
Shape_boundary_Spenalized<-function(x_j,vect_w,j,ul){
  wj<-vect_w[j]
  rest_<-ul-(x_j*vect_w[-j])
  xj<-rest_/wj
  lines(x_j,xj,col="blue")
}
weightFun_sum_penalized<-function(x,u,vect_w){
  return(x * (1 - exp(1-sum_penalized(x=x/u,vect_w = vect_w))))
}
DweightFun_sum_penalized<-function(x,u,vect_w){
  return((1 - exp(1-sum_penalized(x=x/u,vect_w = vect_w)))+ 
           (x*vect_w/u)*exp(1-sum_penalized(x=x/u,vect_w = vect_w)))
}
# Q_norm with weights ------------------------------------------------------------------
########
Norm_q<-function(x,q,weights_nq){
  return((sum(weights_nq*(x)^q)^(1/q)))
}
DNorm_q<-function(x,q,weights_nq){
  Constant<-sum(weights_nq*(x)^q)^((1/q)-1)
  Deriv<-weights_nq*(x)^(q-1)
  return(Constant*Deriv)
}
weightFun_nq<-function(x,u,q,weights_nq){
  N_q<-Norm_q(x/u,q,weights_nq)
  return(x * (1 - exp(-(N_q-1))))
}
DweightFun_nq<-function(x,u,q,weights_nq){
  N_q<-Norm_q(x/u,q,weights_nq)
  First<-(1 - exp(-(N_q - 1))) 
  Second<-(x/u)*DNorm_q(x/u,q,weights_nq)* exp(- (N_q - 1))
  return(First+Second)
}

Wrapper_function<-function(coord,Gam,d,k_taken){
  return(Gamma_to_Sigma_ij(i = coord[1],j = coord[2],
                           Gamma_matrix = Gam,
                           d = d,k_taken=k_taken) )
}
# Obtain the desired input of mgp simulator -------------------------------
### With 0 diagonals and lambda terms.
Gamma_to_Sigma_ij<-function(i,j,Gamma_matrix,d,k_taken){
  orig_coords<-c(1:d)
  first_elt<-Gamma_matrix[i,k_taken]
  second_elt<-Gamma_matrix[j,k_taken]
  third_elt<-Gamma_matrix[i,j]
  #divide by 2 if variogram, by 1 if semi.
  Num<-first_elt+second_elt-third_elt
  return(Num/2)
}
simulation_mixture_HR<-function(d,r,Gamma_mat,A){
  w<-mass_of_scenario(d, r, Sigma, A)  
  T<-rep(-Inf,d)
  b<-sample(c(1:r), prob=w,size=1)
  sign_column_b<-which(A[,b]>0)
  n_column_b<-(A[sign_column_b,b])/sum(A[sign_column_b,b])  
  accept=FALSE
  while(!(accept)){
    if(length(sign_column_b)==1){a<-sign_column_b}
    else{a<-sample(sign_column_b, prob=n_column_b,size=1)}
    sig <- which(A[,b]>0)
    lsig <- length(sig)
    if(lsig == 1){ Z[[k]][sig , ] <- -1/ log( runif(N , 0 , 1) ) * A[sign_column_b, b] 
    } 
    else{
      sub_sigma <- sigma[sig , sig]
    }
    T[sign_column_b]<- t(rmev(N, lsig, sigma = sub_sigma, model = "hr"))*A[sign_column_b, b]
    #T[sign_column_b]<- HR_generator_from_Gamma(a ,b, A[sign_column_b, b], sign_column_b, Sigma[[b]],r)  
    U_0<-runif(1,min=0,max=1)
    if ( (U_0 < rejection_sampling(T)) ) {
      accept = TRUE
    }
  }
  E<-rexp(1,1)
  Y <- T - max(T) + E
}
HR_generator_from_Gamma<-function(i,j, column, sign, covariance_matrix, r){
  # Formula with j= (i of the article).
  # covariance matrix= (sigma of the article)
  T<-mvrnorm(1, mu = -(1 / 2) * diag(covariance_matrix)[sign] + log(r * column) + covariance_matrix[sign, i], 
             Sigma = covariance_matrix[sign, sign])
  
  return(T)
} 