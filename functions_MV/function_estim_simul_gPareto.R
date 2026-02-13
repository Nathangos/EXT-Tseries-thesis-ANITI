
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

Simul_RF_max_NP<-function(Params_risk_Function,M,QTH,
                          Vect_l_function){
  Lg<-apply(X =Vect_l_function,MARGIN = 1,FUN = max)
  Threshold_exp<-quantile(Lg,QTH)
  Inds_exts_lg<-which(Lg>Threshold_exp)
  Excedents_lprime<-Vect_l_function[Inds_exts_lg,]
  # Simulations with LEGRAND method -------------------------------------------------------
  Sim_l<-Non_param_LEGRAND_RiskF(Data_scale_exp =Excedents_lprime,
                                   nb_simul = M,
                                   Threshold_EXP=Threshold_exp)
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
Non_param_LEGRAND_RiskF<-function(Data_scale_exp,nb_simul,Threshold_EXP){
  # From algorithm 1 in LEGRAND ---------------------------------------------
  # Simul Z exp -------------------------------------------------------------
  Zsimul<-rexp(n = nb_simul,rate = 1)
  D<-ncol(Data_scale_exp)
  
  # Simul T  ----------------------------------------------------------------
  Simul_Z<-matrix(Zsimul,ncol = D,
                  nrow=nb_simul,
                  byrow = FALSE)
  print(Data_scale_exp[2,])
  Delta_i<-as.numeric(t(diff(t(Data_scale_exp))))*(-1)
  print(Delta_i[2])
  Max_Vl<-apply(X = Data_scale_exp,MARGIN = 1,FUN = max)
  
  # Independence test BWN ang and rad ---------------------------------------
  print("Test independence")
  print(cor.test(x = Delta_i,y=Max_Vl,method="kendall"))
  print(cor.test(x = Delta_i,y=Max_Vl,method="spearman"))
  plot(Delta_i,Max_Vl,main="Distribution max versus Difference")
  print("Cop Indep test. Null hypothesis = independence")
  UNIF_for_test<-VineCopula::pobs(cbind(Delta_i,Max_Vl))
  IndepCop_chosen<-VineCopula::BiCopIndTest(UNIF_for_test[,1],
                                      UNIF_for_test[,2])
  print(IndepCop_chosen)
  
  # Sample Delta ------------------------------------------------------------
  Delta_tilde<-sample(Delta_i,size = nb_simul,
                      replace = TRUE)
  IND_plus<-which(Delta_tilde>=0)
  Delta_mat<-matrix(0,ncol = D,
                    nrow=nb_simul)
  Delta_mat[-IND_plus,1]<-Delta_tilde[-IND_plus]
  Delta_mat[IND_plus,2]<-(-1)*Delta_tilde[IND_plus]
  
  # Simul complement --------------------------------------------------------
  Simul_RiskF<-Simul_Z+Delta_mat
  return(Simul_RiskF)
  # Exp to orig -------------------------------------------------------------
  # Simul_Z_Unif<-1-exp(-Simul_RiskF)
  # # Simul_RiskF<-Simul_Z_Unif
  # # for(Z in c(1:ncol(Simul_RiskF))){
  # #   Orig_z<-Simul_Z_Unif[,Z]
  # #   Theta_EXTGP<-MOD_EGPD_RiskF[[Z]]
  # #   Orig_to_EGPD<-qEGPDModel1(p = Orig_z,mu = Theta_EXTGP[["mu"]],
  # #                             sigma = Theta_EXTGP[["sigma"]],
  # #                             nu = Theta_EXTGP[["nu"]])
  # #   Simul_RiskF[,Z]<-Orig_to_EGPD
  # # }
  # return(Simul_RiskF)
}

#' Title
#'
#' @param Params_risk_Function 
#' @param M 
#' @param Q_thresh
#' @param Vect_l_function 
#'
#' @return Estimated parameter for the g Pareto law
#' @export
#'
#' @examples
Estim_param_RF_homogeneous<-function(Params_risk_Function,Q_thresh,Vect_l_function,
                                     d){
  Lg<-apply(X =Vect_l_function,MARGIN = 1,
            FUN = Params_risk_Function[["function"]])
  Threshold<-quantile(Lg,Q_thresh)
  Inds_exts_lg<-which(Lg>Threshold)
  Excedents_lprime<-as.matrix(Vect_l_function[Inds_exts_lg,])
  if(Params_risk_Function[["name_RF"]]=="max"){
    l_AIC<-list()
    l_models<-list()
    vect_models<-c("log", "alog", "hr", "neglog", "aneglog", "bilog",
      "negbilog", "ct", "amix")
    for(nameModel in vect_models){
      sub_ob<-evd::fbvpot(x = Vect_l_function,
                          threshold = rep(Threshold,2),
                          model = nameModel)
      l_AIC[[nameModel]]<-AIC(sub_ob)
      l_models[[nameModel]]<-sub_ob
    }
    Model<-Params_risk_Function[["name_model"]]
    BV_method<-evd::fbvpot(x = Vect_l_function,
                        threshold = rep(Threshold,2),
                        model = Model)
    ALL_params<-BV_method$estimate
    Alpha_Likely<-ALL_params[c(5:length(ALL_params))]
    return(list("alpha"=Alpha_Likely,
                "result_AIC"=l_AIC,
                "other"=ALL_params))
    
  }else{
    # Set up, use of g and definition of the weight function  --------------------------
    Model<-Params_risk_Function[["name_model"]]
    Excedents_lprime<-as.list(as.data.frame(t(Excedents_lprime)))
    weightFun<-Params_risk_Function[["weightFunction"]]
    DweightFun<-Params_risk_Function[["dWeightFunction"]]
    
    # Optimisation of the objective function --------------------------------------------
    if(Model=="log"){
      Optimization_Score<-optim(par = Params_risk_Function[["init_opt_param"]],
                                fn = Objective_function,weightFun = weightFun ,
                                lower=0,upper=1,method = "Brent",
                                dWeightFun = DweightFun,Extreme_inds  = Excedents_lprime,
                                u=Threshold,
                                model_used=Model,
                                d=d)
      
    }
    if(Model=="bilog"){
      Optimization_Score<-optim(par =Params_risk_Function[["init_opt_param"]],
                                fn = Objective_function,weightFun = weightFun ,
                                method = "Nelder-Mead",
                                dWeightFun = DweightFun,Extreme_inds  = Excedents_lprime,
                                u=Threshold,
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
                                u=Threshold,model_used=Model,
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
  if(model_used=="bilog"){
    cond_1<-any(theta>=1)
    cond_2<-any(theta<=0)
    if(cond_1||cond_2){
      return(return(1e50))
    }
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
  if(model_used=="bilog"){
    ALPHA<-theta[1]
    BETA<-theta[2]
    R<-sum(PTO_ind)
    W<-PTO_ind[1]/R
    U_Value<-u_bilog_function(alpha = ALPHA,beta = BETA,
                              w = W)
    if((U_Value==0)||(U_Value==1)){
      print(U_Value)
    }
    gradient<-bilog_gradient(alpha =ALPHA,beta = BETA,
                              vect_x = PTO_ind,
                             u_value = U_Value)
    diagHessian<-bilog_hessian_diag(alpha = ALPHA,beta = BETA,
                                       vect_x = PTO_ind,
                              u_value =  U_Value)
    Score_i<-sum(2 * (weights * dWeights) * gradient + weights^2 * 
                 diagHessian + 1/2 * weights^2 * gradient^2)
    return(Score_i)
  }
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
#' Title
#'
#' @param alpha 
#' @param beta 
#' @param vect_x 
#'
#' @return
#' @export
#'
#' @examples
bilog_gradient<-function(alpha,beta,vect_x,u_value){
  ### less doubts. 
  d<-length(vect_x)
  coeff_rad<-(-1)*(d+1)
  r<-sum(vect_x)
  x1<-vect_x[1]
  w<-x1/r
  grad<-rep(NA,length(vect_x))
  dw_x1<-(r-x1)*r^(-2)
  dw_x2<-(-x1)*r^(-2)
  first_member<-coeff_rad*rep((r)^(-1),length(vect_x))
  angle_deriv<-angle_deriv_bilog(alpha,beta,w,
                                 u_w_alpha_beta =u_value)
  second_member<-angle_deriv*c(dw_x1,dw_x2)
  grad<-first_member+second_member
  return(grad)
}
bilog_hessian_diag<-function(alpha,beta,vect_x,u_value){
  ### To check. 
  d<-length(vect_x)
  coeff_rad<-(-1)*(d+1)
  diag_hess<-rep(NA,d)
  x1<-vect_x[1]
  r<-sum(vect_x)
  w<-(x1)*r^(-1)
  deriv_spectral<-angle_deriv_bilog(alpha = alpha,beta=beta,
                                    w =w,u_w_alpha_beta = u_value)
  second_deriv_spectral<-angle_second_deriv_bilog(alpha = alpha,beta=beta,
                                                  w =w,u_value=u_value)
  dw_x1<-(r-x1)*r^(-2)
  dw_x2<-(-x1)*r^(-2)
  dlambda_dr<-coeff_rad*r^(-2)*(-1+deriv_spectral*(w-1))
  dlambda_dw<-r^(-1)*(-deriv_spectral+second_deriv_spectral*(1-w))
  dr_x1<-1
  dr_x2<-1
  ### First coordinate
  first_member<-dlambda_dr*dr_x1
  second_member<-dlambda_dw*dw_x1
  diag_hess[1]<-first_member+second_member
  ### Second coordinate
  dlambda_2_dr<-coeff_rad*r^(-2)*(-1+deriv_spectral*w)
  dlambda_2_dw<--r^(-1)*(deriv_spectral+w*second_deriv_spectral)
  first_member_x2<-dlambda_2_dr*dr_x2
  second_member_member_x2<-dlambda_2_dw*dw_x2
  diag_hess[2]<-first_member_x2+second_member_member_x2
  return(diag_hess)
}

### in the following, deriv of the spectral density

#' Title
#'
#' @param alpha 
#' @param beta 
#' @param w 
#'
#' @return
#' @export
#'
#' @examples
angle_deriv_bilog<-function(alpha,beta,w,u_w_alpha_beta){
  ### Correct formula
  first_part<-(1-w)^(-1)-2*(w)^(-1)
  second_part<-0
  third_part<-0
  if((u_w_alpha_beta!=0)&(u_w_alpha_beta!=1)){
    dv_u<-u_bilog_deriv(alpha = alpha,
                        beta = beta,
                        w = w,u_value = u_w_alpha_beta)
    second_part<-dv_u*(1-u_w_alpha_beta)^(-1)+(1-alpha)*dv_u*(u_w_alpha_beta)^(-1)
    denom_third<-alpha+u_w_alpha_beta*(-alpha+beta)
    third_part<-dv_u*(-alpha+beta)*(denom_third)^(-1)
  }
  First_deriv_spect<-first_part+second_part-third_part
  return(First_deriv_spect)
}

u_bilog_function<-function(alpha,beta,w){
  ### Formula from Coles.
  gmafn <- function(z) (1 - alpha) * (w) * (1 - z)^beta -
    (1 - beta) * (1-w) * z^alpha
  if (is.na(w))
    gma<- NA
  else if (w == 0)
    gma<- 0
  else if (w== 1)
    gma <- 1
  else gma <- uniroot(gmafn, lower = 0, upper = 1, tol = .Machine$double.eps^0.5)$root
  return(gma)
}
u_bilog_function_craft<-function(alpha,beta,w){
  ### Formula from Coles.
  gmafn <- function(z) (1 - alpha) * (1-w) * (1 - z)^beta - 
    (1 - beta) * (w) * z^alpha
  if (is.na(w)) 
    gma<- NA
  else if (w == 0) 
    gma<- 1
  else if (w== 1) 
    gma <- 0
  else gma <- uniroot(gmafn, lower = 0, upper = 1, tol = .Machine$double.eps^0.5)$root
  return(gma)
}
u_bilog_deriv<-function(alpha,beta,w,u_value){
  ### To check !
  if((u_value==0)||(u_value==1)){
    return(0)
  }else{
    first_member<-((1-u_value)*u_value)
    second_member<-(w*(1-w))^(-1)
    third_member<-(u_value*(alpha-beta)-alpha)^(-1)
    return(first_member*second_member*third_member)
  }
  
}
u_bilog_second_deriv<-function(alpha,beta,w,u_value){
  if((u_value==0)||(u_value==1)){
    return(0)
  }else{
    dv_u<-u_bilog_deriv(alpha = alpha,beta=beta,
                        w=w,u_value)
    #(u_value*(-beta+alpha)-alpha)
    funct_alpha_beta<-(u_value*(beta-alpha)+alpha)
    #(beta-alpha)*(dv_u)
    deriv_funct_alpha_beta<-(-beta+alpha)*(dv_u)
    
    ### Define f' used at each step
    dev_1<-((2*w)-1)*(w*(1-w))^(-2)
    dev_2<-dv_u*(1-2*u_value)
    dev_3<-(-1)*deriv_funct_alpha_beta*funct_alpha_beta^(-2)
    
    ### Complete formula
    first_member<-dev_1*(1-u_value)*(u_value)*funct_alpha_beta^(-1)
    second_member<-(w*(1-w))^(-1)*dev_2*funct_alpha_beta^(-1)
    third_member<-(w*(1-w))^(-1)*(1-u_value)*(u_value)*dev_3
    
    Second_deriv<-first_member+second_member+third_member
    return(Second_deriv)
  }
  
}


angle_second_deriv_bilog<-function(w,alpha,beta,u_value){
  dv_u<-u_bilog_deriv(alpha = alpha,beta=beta,
                      w=w,u_value = u_value)
  dv_second_u<-u_bilog_second_deriv(w=w,alpha=alpha,
                                    beta=beta,u_value = u_value)
  ### First member
  first_member<-(1-w)^(-2)+2*w^(-2)
  second_member<-0
  third_member<-0
  fourth_member<-0
  if((u_value!=0)&(u_value!=1)){
    ### Second member
    second_member_one<-(dv_second_u*(1-u_value)+(dv_u)^(2))
    second_member_two<-(1-u_value)^(2)
    second_member<-second_member_one/second_member_two
    
    ### Third member
    third_member_one<-(1-alpha)*((dv_second_u)*u_value-(dv_u)^(2))
    third_member_two<-u_value^2
    third_member<-third_member_one/third_member_two
    
    ### Fourth member
    fourth_member_denom<-(alpha+(u_value)*(-alpha+beta))^(2)
    fourth_member_one_one<-(-alpha+beta)*dv_second_u*(alpha+u_value*(-alpha+beta))
    fourth_member_one_two<-((-alpha+beta)^(2))*(dv_u**(2))
    fourth_member_one<-fourth_member_one_one-fourth_member_one_two
    fourth_member<-fourth_member_one/fourth_member_denom
  }
  Second_deriv_spect<-first_member+second_member+third_member+fourth_member
  return(Second_deriv_spect)
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
Simulation_gParetoP<-function(Params_risk_Function,d,M,
                              theta_opt){
  Nameg<-Params_risk_Function[["name_RF"]]
  name_model<-Params_risk_Function[["name_model"]]
  if(Nameg%in%c("sum","max","min","l2")){
    if(!name_model%in%c("hr","schlater","br")){
      Simul_g<-mev::rgparp(n = M,shape = 1,risk =Nameg ,
                         d = d,
                         model = name_model,
                         param =theta_opt,
                         loc=rep(1,d),
                         scale=rep(1,d))
    }else{
      Sigma_hr_sub<-matrix(0,nrow = d,ncol = d)
      Sigma_hr_sub[2,1]<-theta_opt
      Sigma_hr<-Sigma_hr_sub+t(Sigma_hr_sub)
      Simul_g<-mev::rgparp(n = M,shape = 1,risk =Nameg ,
                         d = d,
                         model = name_model,
                         sigma = Sigma_hr,
                         loc=rep(1,d),
                         scale=rep(1,d))
    }
  }
  else{
    Simul_g<-replicate(M,Rejection_sampling_gPto(Params_risk_Function=Params_risk_Function,
                                        d=d,
                                        theta_opt=theta_opt))
    Simul_g<-matrix(Simul_g,ncol=d,byrow = TRUE)
  }
  return(Simul_g)
}

#' Rejection_sampling_gPto
#'
#' @param Params_risk_Function: list. Parameters of the simulation.  
#' @param d: int. Number of dimensions. 
#' @param theta_opt: float. Estimated estimator of the copula in
#' the chosen model. 
#'
#' @return Simulations of g pareto process for functions not handled by the package
#' @export
#'
#' @examples
Rejection_sampling_gPto<-function(Params_risk_Function,d,
                                  theta_opt){
  
  kept<-FALSE
  g_chosen<-Params_risk_Function[["function"]]
  Model<-Params_risk_Function[["name_model"]]
  while(kept==FALSE){
    unif_sample<-runif(n = 1)
    Copy<-Params_risk_Function
    # U generator
    if(Model!="hr"){
      # candidate<-mev::rmev(n =1,d = d,param = theta_opt,
      #                  model = Model)
      candidate<-mev::rmevspec(n =1,d = d,param = theta_opt,
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
      sim_pto<-extRemes::revd(n = 1,type = "GP",
                              scale=Threshold,
                              shape = 1)
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
#
#
Norm_q<-function(x,q,weights_nq){
  return((sum(weights_nq*(x)^q)^(1/q)))
}
# 
DNorm_q<-function(x,q,weights_nq){
  # Personal
  # Constant<-sum(weights_nq*(x)^q)^((1/q)-1)
  # Deriv<-weights_nq*(x)^(q-1)
  #return(Constant*Deriv)
  # Work from Legrand, Opitz.
  objet<-sum(weights_nq*x**(q-1))*sum(weights_nq*x**q)**((1/q)-1)
  return(objet)
}

weightFun_nq<-function(x,u,q,weights_nq){
  N_q<-Norm_q(x/u,q,weights_nq)
  return(x * (1 - exp(-(N_q-1))))
}

DweightFun_nq<-function(x,u,q,weights_nq){
  N_q<-Norm_q(x/u,q,weights_nq)
  First<-(1 - exp(-(N_q - 1))) 
  Second<-(x)*DNorm_q(x/u,q,weights_nq)*exp(- (N_q - 1))
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