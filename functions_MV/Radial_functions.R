### Mouradhi functions
param_estim <- function(d, r,  grid, lambda, num_col = NULL, start, type = c("SSR_row_HR", "SSR_row_log"), p, w) {
  q <- nrow(grid)
  if (is.null(num_col)) num_col <- r
  l <- d * num_col
  v <- ifelse(type == "SSR_row_log", 1, (d * (d - 1)) / 2)
  
  grid_flat <- as.double(c(t(grid)))  # precompute
  A_length <- d * num_col
  lower_bounds <- rep(0, l + v)
  # interm function
  interm <- function(theta) {
    theta_A <- theta[1:l]
    if (any(theta_A > 1 | theta_A < 0)) return(1e16)
    
    if (type == "SSR_row_log") {
      theta_alpha <- theta[l + 1]
      if (theta_alpha > 1) return(1e16)
      result <- .C(type, as.double(p), as.double(lambda), as.double(theta_A), as.integer(d),
                   as.integer(num_col), as.integer(q), as.double(rep(theta_alpha, num_col)), as.double(w),
                   grid_flat, R = double(1))$R
    } else {
      theta_Gamma <- theta[(l + 1):(l + v)]
      if (any(theta_Gamma < 0)) return(1e16)
      result <- .C(type, as.double(p), as.double(lambda), as.double(theta_A), as.integer(d),
                   as.integer(num_col), as.integer(q), as.double(theta_Gamma), as.double(w),
                   grid_flat, R = double(1))$R
    }
    
    if (!is.finite(result)) return(1e16)
    result
  }
  
  # Initialization
  start_dep <- if (v == 1) runif(1, 0.1, 0.9) else rep(1 , v)
  start_total <- c(start, start_dep)
  
  temp <- optim(start_total, interm, method = "L-BFGS-B", lower = lower_bounds, control = list(maxit = 1000))
  estim_pr <- temp$par
  theta_A_pr <- estim_pr[1:l]
  A_pr <- matrix(theta_A_pr, nrow = d, byrow = TRUE)
  theta_Z_pr <- estim_pr[(l + 1):(l + v)]
  
  # Reuse variables
  A_vec <- as.vector(t(matrix(normalize_group(theta_A_pr, num_col), ncol = num_col, byrow = TRUE)))
  
  # Partial Estimation of A
  interm_2 <- function(theta) {
    if (any(theta > 1 | theta < 0)) return(1e16)
    dep <- if (type == "SSR_row_log") rep(1, num_col) else theta_Z_pr
    
    result <- .C(type, as.double(p), as.double(lambda), as.double(theta), as.integer(d),
                 as.integer(num_col), as.integer(q), as.double(dep), as.double(w),
                 grid_flat, R = double(1))$R
    
    if (!is.finite(result)) return(1e16)
    result
  }
  temp <- optim(theta_A_pr, interm_2, method = "L-BFGS-B", lower = rep(0, l), control = list(maxit = 1000))
  PLS_A <- temp$par
  matrix_A_PLS <- matrix(normalize_group(PLS_A, num_col), ncol = num_col, byrow = TRUE)
  A_vec_fixed <- as.vector(t(matrix_A_PLS))
  
  # Partial Estimation of θ_Z
  interm_3 <- function(theta) {
    if (any(theta < 0) || (type == "SSR_row_log" && any(theta > 1))) return(1e16)
    result <- .C(type, as.double(p), as.double(lambda), A_vec_fixed, as.integer(d),
                 as.integer(num_col), as.integer(q), as.double(theta), as.double(w),
                 grid_flat, R = double(1))$R
    if (!is.finite(result)) return(1e16)
    result
  }
  temp <- optim(theta_Z_pr, interm_3, method = "L-BFGS-B", lower = rep(0, v), control = list(maxit = 1000))
  theta_Z_final <- if (type == "SSR_row_HR") construct_symmetric_matrix(temp$par) else temp$par
  
  return(list(pls_matrix = matrix_A_PLS, pls_dep = theta_Z_final))
}

starting_point <-function(data, nrcol, quant = 0.9){
  N <- nrow(data)
  
  dataP <- apply(data, 2, function(i) N/(N + 0.5 - rank(i)))
  U <- quantile(rowSums(dataP),quant)
  dataU <-dataP[rowSums(dataP)>U,]
  ndata <- t(apply(dataU,1, function(i) i/sum(i)))
  
  kmean <- kmeans(ndata,centers=nrcol,nstart=5)
  startk <- sapply(c(1:nrcol), function(j) kmean$centers[j,]*kmean$size[j])
  resk <- t(apply(startk, 1, function(x) x/sum(x)))
  return(resk)
}
### Simulations from Mouradi
mgpd_simulation_mixture_HR<-function(d,r,Sigma,A){
  w<-mass_of_scenario(d, r, Sigma, A)  
  T<-rep(-Inf,d)
  b<-sample(c(1:r), prob=w,size=1)
  sign_column_b<-which(A[,b]>0)
  n_column_b<-(A[sign_column_b,b])/sum(A[sign_column_b,b])  
  accept=FALSE
  while(!(accept)){
    if(length(sign_column_b)==1){a<-sign_column_b}
    else{a<-sample(sign_column_b, prob=n_column_b,size=1)}
    T[sign_column_b]<- HR_generator(a ,b, A[sign_column_b, b], 
                                    sign_column_b, Sigma[[b]],r)  
    U_0<-runif(1,min=0,max=1)
    if ( (U_0 < rejection_sampling(T)) ) {
      accept = TRUE
    }
  }
  E<-rexp(1,1)
  Y <- T - max(T) + E
  return(Y)
}

mgpd_simulation_mixture_logistic<-function(d,r,alpha,A, type = c("logistic", "HR")){
  w<-colSums((A)^(1/alpha))^(alpha)/sum(colSums((A)^(1/alpha))^(alpha))   
  T<-rep(-Inf,d)
  b<-sample(c(1:r), prob=w,size=1)
  sign_column_b<-which(A[,b]>0)
  n_column_b<-(A[sign_column_b,b])/sum(A[sign_column_b,b])  
  accept=FALSE
  while(!(accept)){
    if(length(sign_column_b)==1){a<-sign_column_b}
    else{a<-sample(sign_column_b, prob=n_column_b,size=1)}
    T[sign_column_b]<-transformed_logsitic_generator(a, sign_column_b, A[sign_column_b,b], alpha[b],r)
    U_0<-runif(1,min=0,max=1)
    if ( (U_0 < rejection_sampling(T)) ) {accept = TRUE}
  }
  E<-rexp(1,1)
  Y <- T - max(T) + E
  return(Y)
}
Sample_cond_g_Rjection_Sampling<-function(Mu_vector,
                                          Cov_mat,g,Th_g){
  verified<-FALSE
  d<-ncol(Cov_mat)
  while(verified==FALSE){
    vect_norm<-c(mvtnorm::rmvnorm(n = 1,mean = Mu_vector,
                                sigma = Cov_mat))
    g_sim<-g(vect_norm)
    if(g_sim>Th_g){
      verified<-TRUE
    }
  }
  return(vect_norm)
  
}
Sample_cond_Laplacian_Rjection_Sampling<-function(Mu_ALD,
                                          Scatter_ALD,g,Th_g){
  verified<-FALSE
  d<-ncol(Scatter_ALD)
  while(verified==FALSE){
    vect_ald<-c(L1pack::rmLaplace(n = 1,center = Mu_ALD,
          Scatter = Scatter_ALD))
    g_sim<-g(vect_ald)
    if(g_sim>Th_g){
      verified<-TRUE
    }
  }
  return(vect_ald)
  
}

### Simulation from exceedances of Gaussian 
Sample_cond_g<-function(Mu_vector,Cov_mat,
                        g,Th_g){
  verified<-FALSE
  d<-ncol(Cov_mat)
  while(verified==FALSE){
    j_chosen<-sample(c(1:d),size = 1)
    vect_norm<-rep(NA,d)
    W_sim<-truncnorm::rtruncnorm(1,a=Th_g,b=Inf)
    Sigma_inv_W<-solve(Cov_mat[j_chosen,j_chosen])
    Cov_U_condW<-Cov_mat[-j_chosen,j_chosen]
    Rest_mean<-Cov_U_condW%*%Sigma_inv_W%*%(W_sim-Mu_vector[j_chosen])
    mu_vector_cond<-Mu_vector[-j_chosen]+Rest_mean
    Rest_cov_mat<-Cov_U_condW%*%Sigma_inv_W%*%t(Cov_U_condW)
    cov_mat_cond<-Cov_mat[j_chosen,j_chosen]-Rest_cov_mat
    if(d>2){
      
      Simul_rest<-mvtnorm::rmvnorm(
        n = 1,mean = mu_vector_cond,
        sigma =cov_mat_cond)
    }else{
      cov_mat_cond<-as.numeric(cov_mat_cond)
      mu_vector_cond<-as.numeric(mu_vector_cond)
      sig<-sqrt(cov_mat_cond)
      Simul_rest<-rnorm(1,mean = mu_vector_cond,
                        sd = sig)
    }
    vect_norm[-j_chosen]<-Simul_rest
    vect_norm[j_chosen]<-W_sim
    g_sim<-g(vect_norm)
    if(g_sim>Th_g){
      verified<-TRUE
    }
  }
  return(vect_norm)
}
#' Simul_from_Htawn
#'
#' @param model_mex_all: MexALL object.
#' @param n_sim: int. Number of simulations.
#' @param d: int. Number of dimensions.
#' @param Name_vars: list(str)
#' @param prop_class: float. Scenario probability
#' @param ind_ref: int. First case chosen. 
#' @param thresh_censor: float. Threshold used for (-ind_ref) cases.
#' @param Q_sim: float. Quantile used for the conditioning variable.
#'
#' @return
#' @export
#'
#' @examples
Simul_from_Htawn<-function(model_mex_all,n_sim,d,Name_vars,
                           prop_class,ind_ref,
                           thresh_censor,Q_sim){
  Data_rad<-data.frame(matrix(NA,nrow = n_sim,ncol=d))
  colnames(Data_rad)<-Name_vars
  ### same value for prob above threshold.
  list_sample<-list()
  probs_class<-prop_class/sum(prop_class)
  sample_marg<-sample(c(1:d),size = n_sim,replace = TRUE,
                      prob = probs_class)
  for(z in c(1:d)){
    Ind_z<-which(sample_marg==z)
    L_z<-length(Ind_z)
    Model_z<-model_mex_all[[z]]
    Qsim_z<-Q_sim[z]
    n_sim_z<-length(Ind_z)
    Pop_sims_z<-t(replicate(n_sim_z,
                  HTAWN_censorshisp(model_z = Model_z,
                    prop_z = Qsim_z,thresh_censor = thresh_censor,
                    ind_ref = ind_ref,z = z)))
    for(variable in Name_vars){
      Data_rad[Ind_z,variable]<-Pop_sims_z[,variable]
    }
    
  }
  return(Data_rad)
}
HTAWN_censorshisp<-function(model_z,prop_z,thresh_censor,ind_ref,z){
    accept<-FALSE
    while(accept==FALSE){
      Sims_Cond_z<-predict(object= model_z,
                                    pqu=prop_z,nsim =2,
                                    which = z)$data$simulated
      Sims_Cond_z<-unlist(Sims_Cond_z[1,])
      u_z<-Sims_Cond_z[z]
      if(z==ind_ref){
        return(Sims_Cond_z)
      }else{
        if(u_z<thresh_censor){
          return(Sims_Cond_z)
        }
      }
      
    }
}

### HGD + Laplace distributions
#' Convert_HGD_Unif
#'
#' @param N_test: int. Number of simulated vectors.
#' @param Model_HGD: Class. HGD model fitted.
#'
#' @return Vectors from the HGD model with uniform margins
#' @export
#'
#' @examples
Simulations_fromHGD_Copula<-function(N_test,
                                   Model_HGD){
  
  Sim_HGD_scale<-ghyp::rghyp(n = N_test,
                             object=Model_HGD)
  Sim_Unif_conv<-sapply(X = c(1:ncol(Sim_HGD_scale)),
                        function(z,Vectors){
                          sigZ<-Model_HGD@sigma[z,z]
                          univariate.ghyp <- ghyp::ghyp(
                            lambda = Model_HGD@lambda,  
                            psi = Model_HGD@psi,
                            chi = Model_HGD@chi,
                            sigma = sigZ)
                          if(is.null(dim(Vectors))){
                            Obs_z<-Vectors[z]
                          }else{
                            Obs_z<-Vectors[,z]
                          }
                          return(ghyp::pghyp(Obs_z,
                                             object = univariate.ghyp))
                        },
                        Vectors=Sim_HGD_scale)
  return(Sim_Unif_conv)
  
}
Sample_cond_HGD_Rjection_Sampling<-function(Model_HGD,g,Th_g){
  verified<-FALSE
  d<-ncol(Model_HGD@sigma)
  while(verified==FALSE){
    Copula_HGD<-Simulations_fromHGD_Copula(N_test = 1,
                Model_HGD = Model_HGD)
    vect_N<-qnorm(Copula_HGD)
    g_sim<-g(vect_N)
    if(g_sim>Th_g){
      verified<-TRUE
    }
  }
  return(vect_N)
  
}

Fit_HGD_from_Laplace<-function(Lap_vectors,Model_Lap_inits){
  
  d<-ncol(Lap_vectors)
  ### Use Generalized Laplace
  Mod_HGD<-ghyp::fit.VGmv(data = Lap_vectors,
                         gamma=rep(0,d),
                         lambda = 1,opt.pars = c(lambda=TRUE,
                                                 mu=TRUE,gamma=FALSE,
                                                 sigma=TRUE),
                         symmetric = TRUE,
                         sigma=Model_Lap_inits$Scatter,
                         mu=Model_Lap_inits$center,
                         standardize=TRUE)
  return(Mod_HGD)
  
}
Convert_marg_model_ghyp<-function(vect,model_univ){
  return(ghyp::qghyp(vect,object = model_univ))
}
