# Approach the cdf of marginals  ------------------------------------------
# knowing the maximum is extreme ------------------------------------------

# Functions for the couple of l -------------------------------------------

#' Function_lower_distrib
#'
#' @param x : float. Quantile of the law.
#' @param ASY : vector(float). Asymmetric parameters. 
#' @param Th : float. Threshold of extreme individuals. 
#' @param DEP : float. Dependence parameter. 
#' @param prob_given : float (between 0 and 1). Sampled value from an uniform.
#'
#' @return Difference between the cdf and a sampled value of an uniform.
#' @export
#'
#' @examples
Function_lower_distrib<-function(x,ASY,Th,DEP,prob_given,z,Model,Cste){
  vect_g<-rep(Th,2)
  vect_g[z]<-x
  # Num<-exp(-(x/Cste[z])^(-1))-evd::pbvevd(vect_g,
  #                                         model=Model,
  #                                         mar1 = c(Cste[1],Cste[1],1),
  #                                         mar2=c(Cste[2],Cste[2],1),
  #                                         dep=DEP,asy=ASY)
  vect_end<-rep(Inf,2)
  vect_end[z]<-x
  Num_a<-1+log(evd::pbvevd(vect_end,
                   model=Model,
                   mar1 = c(Cste[1],Cste[1],1),
                   mar2=c(Cste[2],Cste[2],1),
                   dep=DEP,asy=ASY))
  Num_b<-1+log(evd::pbvevd(vect_g,
                          model=Model,
                          mar1 = c(Cste[1],Cste[1],1),
                          mar2=c(Cste[2],Cste[2],1),
                          dep=DEP,asy=ASY))
  Num<-Num_a-Num_b
  #1-
  Denum<--log(evd::pbvevd(c(Th,Th),
                       mar1 = c(Cste[1],Cste[1],1),
                       mar2=c(Cste[2],Cste[2],1)
                       ,model=Model,
                       dep=DEP,asy=ASY))
  return((Num/Denum)-prob_given)
}


#' Simulate_law_cond_max
#' 
#' @param NB_simulations : int. Desired number of simulations. 
#' @param Params_Biv : list. Parameters of the joint law and constants. 
#' @param u : float. Realisation of the uniform law. 
#' @param Th : float. Threshold outreached by the maximum. 
#' @param z : int. Index of the marginal used. 
#'
#' @return vector(float). (NB_simulations) realisations of the marginal law.
#' @export
#'
#' @examples
Simulate_law_cond_max<-function(Nb_simulations,Params_Biv,Th,z){
  
  U<-runif(n = Nb_simulations)
  S<-sapply(U,FUN = Function_mixture,Params_Biv=Params_Biv,
            z=z,Th=Th)
  return(list("values_F"=S,
              "Unif"=U))
}

#' Function_mixture
#'
#' @param Params_Biv : list. Parameters of the joint law and constants. 
#' @param u : float. Realisation of the uniform law. 
#' @param Th : float. Threshold outreached by the maximum. 
#' @param z : int. Index of the marginal used. 
#'
#' @return Float. Realisation of the marginal law. 
#' @export
#'
#' @examples
Function_mixture<-function(Params_Biv,u,Th,z){
  
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Cste<-Params_Biv$cste
  Model<-Params_Biv$model
  # NUM_QU<-exp(-(Th/Cste[z])^(-1))-evd::pbvevd(c((Th),(Th)),
  #                                             mar1 = c(Cste[1],Cste[1],1),
  #                                             mar2=c(Cste[2],Cste[2],1),
  #                                             model=Model,
  #                                             dep=dep,asy=ASY)
  vect_inf<-rep(Inf,2)
  vect_inf[z]<-Th
  Num_QUa<-1+log(evd::pbvevd(vect_inf,model=Model,
                             mar1 = c(Cste[1],Cste[1],1),
                             mar2=c(Cste[2],Cste[2],1),
                             dep=dep,asy=ASY))
  Num_QUb<-1+log(evd::pbvevd(rep(Th,2),model=Model,
                            mar1 = c(Cste[1],Cste[1],1),
                            mar2=c(Cste[2],Cste[2],1),
                            dep=dep,asy=ASY))
  NUM_QU<-Num_QUa-Num_QUb
  
  # without Taylor ----------------------------------------------------------
  # DENOM_QU<-1-evd::pbvevd(c((Th),(Th)),model=Model,
  #                         mar1 = c(Cste[1],Cste[1],1),
  #                         mar2=c(Cste[2],Cste[2],1),
  #                         dep=dep,asy=ASY)
  DENOM_QU<--log(evd::pbvevd(c((Th),(Th)),model=Model,
                             mar1 = c(Cste[1],Cste[1],1),
                             mar2=c(Cste[2],Cste[2],1),
                             dep=dep,asy=ASY))
  QU<-NUM_QU/DENOM_QU
  if(u<QU){
    
    # Uniroot to find a correct sampled value ---------------------------------
    # knowing the realisation of the uniform ----------------------------------
    
    sampled_value<-uniroot(f = Function_lower_distrib,interval = c(0,Th),
                           ASY=ASY,Th=Th,DEP=dep,prob_given=u,z=z,
                           Cste=Cste,Model=Model)$root
  }
  else{
    
    # Use the upper side of the density ---------------------------------------
    # v<-u*DENOM_QU+evd::pbvevd(c((Th),(Th)),
    #                           mar1 = c(Cste[1],Cste[1],1),
    #                           mar2=c(Cste[2],Cste[2],1),
    #                           model=Model,
    #                           dep=dep,asy=ASY)
    v<-u*DENOM_QU+log(evd::pbvevd(c((Th),(Th)),
                                  mar1 = c(Cste[1],Cste[1],1),
                                  mar2=c(Cste[2],Cste[2],1),
                                  model=Model,
                                  dep=dep,asy=ASY))
    
    #sampled_value<-(-Cste[z])/log(v)
    sampled_value<-uniroot(f = Theor_counterpart_marg,interval = c(Th,10^10),
                           Params_Biv=Params_Biv,Th=Th,z=z,
                           prob_given=v)$root
  }
  return(sampled_value)
}

Theor_counterpart_marg<-function(Params_Biv,Th,z,x,prob_given){
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Cste<-Params_Biv$cste
  Model<-Params_Biv$model
  
  vect_end_marg<-rep(Inf,2)
  vect_end_marg[z]<-x
  return(log(evd::pbvevd(vect_end_marg,
                               mar1 = c(Cste[1],Cste[1],1),
                               mar2=c(Cste[2],Cste[2],1),
                               model=Model,
                               dep=dep,asy=ASY))-prob_given)
}
# Law to uniform ----------------------------------------------------------


CDF_1d_Mbiv<-function(value,Th,Params_Biv,z){
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Cste<-Params_Biv$cste
  Model<-Params_Biv$model
  fct_replace_min<-function(x_value,ref_value){
    return(min(x_value,ref_value))
  }
  if(Model%in%c("log","bilog","hr")){
    Denom<--log(evd::pbvevd(q = rep(Th,2),
                            mar1 = c(Cste[1],Cste[1],1),
                            mar2=c(Cste[2],Cste[2],1),
                            model=Model,
                            alpha=dep[1],beta=dep[2],
                            dep=dep[1],
                            asy=ASY))
    vector_inf<-rep(Inf,2)
    vector_inf[z]<-value
    vector_th<-rep(Th,2)
    vector_th[z]<-value
    if(value<Th){
      First_elt<-log(evd::pbvevd(q=vector_inf,
                                 mar1=c(Cste[1],Cste[1],1),
                     mar2=c(Cste[2],Cste[2],1),
                     model=Model,
                     alpha=dep[1],beta=dep[2],
                     dep=dep[1],
                     asy=ASY))
      
      Second_elt<--log(evd::pbvevd(q=vector_th,
                                 mar1=c(Cste[1],Cste[1],1),
                     mar2=c(Cste[2],Cste[2],1),
                     model=Model,
                     alpha=dep[1],beta=dep[2],
                     dep=dep[1],
                     asy=ASY))
      # Third_elt<--log(evd::pbvevd(q=c(Bound_inf,Inf),
      #                              mar1=c(Cste[1],Cste[1],1),
      #                  mar2=c(Cste[1],Cste[1],1),
      #                  model=Model,
      #                  alpha=dep[1],beta=dep[2],
      #                  asy=ASY))
      Third_elt<-0
      # Fourth_elt<-log(evd::pbvevd(q=c(Bound_inf,Th),
      #                             mar1=c(Cste[1],Cste[1],1),
      #                 mar2=c(Cste[1],Cste[1],1),
      #                 model=Model,
      #                 alpha=dep[1],beta=dep[2],
      #                 asy=ASY))
      Fourth_elt<-0
      Numerator<-First_elt+Second_elt+Third_elt+Fourth_elt
      return(Numerator/Denom)
    }
    else{
      First_elt<--log(evd::pbvevd(q=rep(Th,2),
                                  mar1=c(Cste[1],Cste[1],1),
                                  mar2=c(Cste[2],Cste[2],1),
                                  model=Model,
                                  alpha=dep[1],beta=dep[2],
                                  dep=dep[1],
                                  asy=ASY))
      # Second_elt<--log(evd::pbvevd(q=c(Bound_inf,Inf),
      #                              mar1=c(Cste[1],Cste[1],1),
      #                  mar2=c(Cste[1],Cste[1],1),
      #                  model=Model,
      #                  alpha=dep[1],beta=dep[2],
      #                  asy=ASY))
      Second_elt<-0
      # Third_elt<-log(evd::pbvevd(q=c(Bound_inf,Th),
      #                             mar1=c(Cste[1],Cste[1],1),
      #                 mar2=c(Cste[1],Cste[1],1),
      #                 model=Model,
      #                 alpha=dep[1],beta=dep[2],
      #                 asy=ASY))
      Third_elt<-0
      Fourth_elt<-log(evd::pbvevd(q=vector_inf,
                                  mar1=c(Cste[1],Cste[1],1),
                      mar2=c(Cste[2],Cste[2],1),
                      model=Model,
                      alpha=dep[1],beta=dep[2],
                      dep=dep[1],
                      asy=ASY))
      # Fifth_elt<--log(evd::pbvevd(q=c(x,Bound_inf),
      #                            mar1=c(Cste[1],Cste[1],1),
      #                mar2=c(Cste[1],Cste[1],1),
      #                model=Model,
      #                alpha=dep[1],beta=dep[2],
      #                asy=ASY))
      Fifth_elt<-0
      # Sixth_elt<-log(evd::pbvevd(q=c(Th,Bound_inf),
      #                            mar1=c(Cste[1],Cste[1],1),
      #                mar2=c(Cste[1],Cste[1],1),
      #                model=Model,
      #                alpha=dep[1],beta=dep[2],
      #                asy=ASY))
      Sixth_elt<-0
      Numerator<-First_elt+Second_elt+Third_elt+
        Fourth_elt+Fifth_elt+Sixth_elt
      return(Numerator/Denom)
    }
    # First_elt<-log(evd::pbvevd(q = vector,
    #                            mar1 = c(Cste[1],Cste[1],1),
    #                            mar2=c(Cste[2],Cste[2],1),
    #                            model=Model,
    #                            alpha=dep[1],beta=dep[2],
    #                            asy=ASY))
    # Second_elt<-log(evd::pbvevd(q = Min_,
    #                             mar1 = c(Cste[1],Cste[1],1),
    #                             mar2=c(Cste[2],Cste[2],1),
    #                             model=Model,
    #                             alpha=dep[1],beta=dep[2],
    #                             asy=ASY))
  }
  else{
    # First_elt<-log(evd::pbvevd(q = vector,
    #                            mar1 = c(Cste[1],Cste[1],1),
    #                            mar2=c(Cste[2],Cste[2],1),
    #                            model=Model,
    #                            dep=dep,asy=ASY))
    # Second_elt<-log(evd::pbvevd(q = Min_,
    #                             mar1 = c(Cste[1],Cste[1],1),
    #                             mar2=c(Cste[2],Cste[2],1),
    #                             model=Model,
    #                             dep=dep,
    #                             asy=ASY))
    Denom<--log(evd::pbvevd(q = c(Th,2),
                            mar1 = c(Cste[1],Cste[1],1),
                            mar2=c(Cste[2],Cste[2],1),
                            model=Model,
                            dep=dep,
                            asy=ASY))
  }
  return((First_elt-Second_elt)/Denom)
}


CDF_law_2d<-function(Th,Params_Biv,vect_xy){
  return(CDF_1d_Mbiv(x = vect_xy[1],
                     y = vect_xy[2],
                     Th = Th,
                     Params_Biv = Params_Biv))
  
}

CDF_law_2d_deriv_z<-function(Th,Params_Biv,vect_xy,z){
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Cste<-Params_Biv$cste
  Model<-Params_Biv$model
  Prob_inf<-1+log(evd::pbvevd(q = c(Th,Th),
                        mar1 = c(Cste[1],Cste[1],1),
                        mar2=c(Cste[2],Cste[2],1),
                        model=Model,
                        dep=dep,asy=ASY))
  #-1
  Prob_sup<-1-Prob_inf
  
  # first case--one couple changed, one element becomes null -------------------------------------------------------------
  x_dom<-vect_xy[1]>Th
  if(x_dom){
    j_dom<-1
    XInput<-vect_xy
    XInput[j_dom]<-Th
    #evd::pbvevd(q = vect_xy,
    # mar1 = c(Cste[1],Cste[1],1),
    # mar2=c(Cste[2],Cste[2],1),
    # model=Model,
    # dep=dep,asy=ASY)*
    First<-V_x_y_deriv_z_ind(vect_x_y = vect_xy, z = z,
                              Cst_ = Th,Params_Biv = Params_Biv)
    #evd::pbvevd(q = XInput,
    # mar1 = c(Cste[1],Cste[1],1),
    # mar2=c(Cste[2],Cste[2],1),
    # model=Model,
    # dep=dep,asy=ASY)
    Second<-V_x_y_deriv_z_ind(vect_x_y = XInput,z = z,
                               Cst_ = Th,Params_Biv = Params_Biv)
    Partial_1<-First-Second
  }
  
  # second case --------------------------------------------------------------
  y_dom<-vect_xy[2]>Th
  if(y_dom){
    j_dom<-2
    XInput<-vect_xy
    XInput[j_dom]<-Th
    # evd::pbvevd(q = vect_xy,
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)*
    First<-V_x_y_deriv_z_ind(vect_x_y = vect_xy, z = z,
                              Cst_ = Th,
                              Params_Biv = Params_Biv)
    # evd::pbvevd(q = XInput,
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)*
    Second<-V_x_y_deriv_z_ind(vect_x_y = XInput, z = z,
                               Cst_ = Th,
                               Params_Biv = Params_Biv)
    Partial_2<-(First-Second)
  }
  
  # third case --------------------------------------------------------------
  Min<-min(vect_xy)
  min_dom<-Min>Th
  if(min_dom){

    # without log transf ------------------------------------------------------

    # evd::pbvevd(q = vect_xy,
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)* 
    First<-V_x_y_deriv_z_ind(vect_x_y = vect_xy, z = z,
                            Cst_ = Th,Params_Biv = Params_Biv)
    # evd::pbvevd(q = c(vect_xy[1],Th),
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)*
    Second<-V_x_y_deriv_z_ind(vect_x_y = c(vect_xy[1],Th),
                                                           z = z,
                                                           Cst_ = Th,
                                                           Params_Biv = Params_Biv)

    # comment -----------------------------------------------------------------

    # l -----------------------------------------------------------------------

    # evd::pbvevd(q = c(Th,vect_xy[2]),
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)*
    Third<-V_x_y_deriv_z_ind(vect_x_y = c(Th,vect_xy[2]),
                                                          z = z,
                                                          Cst_ = Th,
                                                          Params_Biv = Params_Biv)
    
    # evd::pbvevd(q =rep(Th,length(vect_xy)),
    #             mar1 = c(Cste[1],Cste[1],1),
    #             mar2=c(Cste[2],Cste[2],1),
    #             model=Model,
    #             dep=dep,asy=ASY)*
    Fourth<-V_x_y_deriv_z_ind(vect_x_y = rep(Th,length(vect_xy)),
                                                           z = z,
                                                           Cst_ = Th,
                                                           Params_Biv = Params_Biv)
    Whole<-First-(Second+Third)+Fourth
  }
  Numtor<-0
  if(x_dom){
    Numtor<-Partial_1
  }
  if(y_dom){
    Numtor<-Partial_2
  }
  if(min_dom){
    Numtor<-Whole
  }
  return(Numtor/Prob_sup)
}


Density_deriv_z<-function(Th,Params_Biv,variable,z){
  
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Cste<-Params_Biv$cste
  Model<-Params_Biv$model
  
  cas_1<-(variable>Th)
  Min_th<-min(variable,Th)
  
  # Inf for considering only one marginal -----------------------------------
  Normaltion_constant<-Params_Biv$cste[z]
  #exp(-(Min_th/Cste[z])^(-1))*
  vect_inf_marg<-rep(Inf,2)
  vect_inf_marg[z]<-Min_th
  # Whole_2<-V_x_deriv_z_ind(x = Min_th,Normaltion_constant=Normaltion_constant,
  #                          Cst_=Th)
  Whole_2<-V_x_y_deriv_z_ind(vect_x_y =  vect_inf_marg,
                             z = z,
                             Cst_ = Th,
                             Params_Biv = Params_Biv)
  vect_given<-rep(Th,2)
  vect_given[z]<-Min_th
  # origin<-evd::pbvevd(q = vect_given,
  #                     mar1 = c(Cste[1],Cste[1],1),
  #                     mar2=c(Cste[2],Cste[2],1),
  #                     model=Model,
  #                     dep=dep,asy=ASY)
  #origin*
  output<-V_x_y_deriv_z_ind(vect_x_y = vect_given,
                                   z = z,
                                   Cst_ = Th,
                                   Params_Biv = Params_Biv)
  Partial_2<-output
  Partial_1<-Whole_2-Partial_2

  # without the log transf --------------------------------------------------
  # comment : ---------------------------------------------------------------

  # exp(-(variable/Cste[z])^(-1)) -------------------------------------------
  
  # vals_x<-V_x_deriv_z_ind(variable,Normaltion_constant=Normaltion_constant,
  #                                                       Cst_=Th)
  vect_inf2<-rep(Inf,2)
  vect_inf2[z]<-variable
  vals_x<-V_x_y_deriv_z_ind(vect_x_y =vect_inf2,
                            z = z,
                            Cst_ = Th,
                            Params_Biv = Params_Biv)
    
  # denominator -------------------------------------------------------------
  
  Prob_inf<-1+log(evd::pbvevd(q = c(Th,Th),
                        mar1 = c(Cste[1],Cste[1],1),
                        mar2=c(Cste[2],Cste[2],1),
                        model=Model,
                        dep=dep,asy=ASY))
  Prob_sup<-1-Prob_inf
  if(cas_1==TRUE){
    Numerator<-vals_x
  }
  else{
    Numerator<-Partial_1
  }
  return(Numerator/Prob_sup)
}


V_x_deriv_z_ind<-function(x,Normaltion_constant,Cst_){
  if(x==Cst_){return(0)}
  else{
    return(Normaltion_constant*x^(-2))
  }
}

V_x_y_deriv_z_ind<-function(vect_x_y,z,Cst_,Params_Biv){
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  z_alt<-ifelse(z==1,2,1)
  Normaltion_constant<-Params_Biv$cste
  Model<-Params_Biv$model
  x<-vect_x_y[1]
  y<-vect_x_y[2]
  
  # Deriv=0 if at index z we have the Cst_ ----------------------------------
  not_null<-!(vect_x_y[z]==Cst_)
  if(not_null==FALSE){
    return(0)
  }
  else{
    if(Model=="alog"){
      Const<-((x/(ASY[1]*Normaltion_constant[1]))^(-1/dep)+
                (y/(ASY[2]*Normaltion_constant[2]))^(-1/dep))^(dep-1)
      Coeff_deriv<-vect_x_y[z]^((-1/dep)-1)*(ASY[z]*Normaltion_constant[z])^(1/dep)
      Coeff_deriv1<-Coeff_deriv*Const
      Value<-Coeff_deriv1+(1-ASY[z])*(vect_x_y[z])^(-2)*Normaltion_constant[z]
    }
    if(Model=="log"){
      Const<-((x/(Normaltion_constant[1]))^(-1/dep)+
                (y/(Normaltion_constant[2]))^(-1/dep))^(dep-1)
      Coeff_deriv<-vect_x_y[z]^((-1/dep)-1)*(Normaltion_constant[z])^(1/dep)
      Value<-Coeff_deriv*Const
    }
    if(Model=="hr"){
      
      ratio<-(vect_x_y[z]/Normaltion_constant[z])/(vect_x_y[z_alt]/Normaltion_constant[z_alt])
      Q_cst<-dep^(-1)+(1/2)*dep*log(ratio^(-1))
      Q_cst2<-dep^(-1)+(1/2)*dep*log(ratio)
      first<-vect_x_y[z]^(-2)*pnorm(q = Q_cst)*Normaltion_constant[z]+
        vect_x_y[z]^(-2)*dnorm(Q_cst)*(1/2)*dep*Normaltion_constant[z]
      second<--vect_x_y[z_alt]^(-1)*dnorm(Q_cst2)*(1/2)*dep*vect_x_y[z]^(-1)
      Value<-first+second
    }
    if(Model=="aneglog"){
      first<-Normaltion_constant[z]*vect_x_y[z]^(-2)
      Cste_<-((vect_x_y[z]/(Normaltion_constant[z]*ASY[z]))^dep
               +(vect_x_y[z_alt]/(Normaltion_constant[z_alt]*ASY[z_alt]))^dep)^((-1/dep)-1)
      second<-Cste_*(-1)*((Normaltion_constant[z]*ASY[z]))^(-dep)*vect_x_y[z]^(dep-1)
      Value<-first+second
    }
    return(Value)
  }
}

Weight_Th<-function(vect_xy,other_variable,z){
  ASY<-Params_Biv$asy
  dep<-Params_Biv$dep
  Normaltion_constant<-Params_Biv$cste
  Model<-Params_Biv$model
  
  x<-vect_xy[1]
  y<-vect_xy[2]
  Const<-((x/(ASY[1]*Normaltion_constant[1]))^(-1/dep)+
            (y/(ASY[2]*Normaltion_constant[2]))^(-1/dep))^(dep-1)
  Coeff_deriv<-vect_x_y[z]^((-1/dep)-1)*(ASY[z]*Normaltion_constant[z])^(1/dep)
  Coeff_deriv1<-Coeff_deriv*Const
  Coeff_cplmt<-Coeff_deriv1+(1-ASY[z])*(vect_x_y[z])^(-2)
  return(Coeff_cplmt)
  
}
Compar_probs<-function(Params_Biv,Th,z,xsimul,prob_given,y){
  Prob_found<-Cond_law_1d(x= xsimul,y = y,Th = Th,Params_Biv = Params_Biv,
                          z = z)
  return(Prob_found-prob_given)
}

Simul_cond_one<-function(input_x,Params_biv,Th,z){
  Cste<-Params_biv$cste
  ASY<-Params_biv$asy
  Model<-Params_biv$model
  dep<-Params_biv$dep
  
  x<-runif(1)
  print(x)
  z_alt<-ifelse(z==1,2,1)
  vect_weight<-rep(Th,2)
  vect_weight[z]<-input_x
  if((input_x<Th)){
    return(uniroot(Compar_probs,interval = c(as.numeric(Th),10^10),xsimul=input_x,
                   Th=Th,Params_Biv=P_BIV,prob_given=x,z=1)$root)
  }
  if((input_x>Th)){
    return(uniroot(Compar_probs,interval = c(0,10^10),xsimul=input_x,
                   Th=Th,Params_Biv=P_BIV,prob_given=x,z=1)$root)
  }
}
