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
  Data_rad<-data.frame(matrix(NA,nrow = n_sim,
                    ncol=d))
  colnames(Data_rad)<-c(Name_vars)
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
        ind_ref = ind_ref,z = z,
        l_name=Name_vars)))
    #return(Pop_sims_z)
    for(variable in Name_vars){
      Data_rad[Ind_z,variable]<-Pop_sims_z[,variable]
    }
    
  }
  return(Data_rad)
}
### Convergence of the defined "angles" of the radial vector
Extreme_covRad_per_K<-function(Df_Rad_Orig,l_name,k,
                                   vect_lg){
  Order_lg<-order(vect_lg,
                  decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  Inds_exts<-which(vect_lg>=Threshold_lg)
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Sub_Rad<-Df_Rad_Orig[Inds_exts,]
  Sub_Lg<-vect_lg[Inds_exts]
  ###
  Theta_of_rad<-t(t(Sub_Rad)%*%diag(Sub_Lg^(-1)))
  Theta_of_rad<-as.data.frame(Theta_of_rad)
  colnames(Theta_of_rad)<-l_name
  Mat_corr<-Fct_correlations(method_corr = "kendall",
                    df1 = Theta_of_rad,
                   Intersect_times = l_name,
                   df2 =NA,mat_corr=FALSE)
  Mat_corr$number_exceed<-rep(k,nrow(Mat_corr))
  return(Mat_corr)
}
Extreme_covRad_evol<-function(Df_Rad_Orig,l_name,vector_k,
                           Ref_RiskF,CPU_hearts){
  All_results<-lapply(X = vector_k,
                 Extreme_covRad_per_K,
                 l_name=l_name,
                 Df_Rad_Orig=Df_Rad_Orig,
                 vect_lg=Ref_RiskF)
  return(do.call(what = rbind.data.frame,
                 All_results))
}
Loop_TexmexCarlo_Extrapol<-function(nSample, mexList, pqu_extrapol,mult = 10){
  d<-length(mexList)
  Matrix_full<-matrix(NA,nrow = nSample,ncol = d)
  N_obtained<-0
  while(N_obtained<nSample){
    Turn_j<-TexmexMCarlo_Extrapol(nSample = nSample,mexList = mexList,
                                  mult = mult,pqu_extrapol = pqu_extrapol)
    Df_j<-do.call("rbind",
                  Turn_j$list_craft)
    CNAMES<-colnames(Df_j)
    Nturn<-nrow(Df_j)
    N_taken<-min(Nturn,nSample)
    
    Gap<-min(N_taken,nSample-N_obtained)
    Beg<-1+N_obtained
    End<-Gap+N_obtained
    Matrix_full[Beg:End,]<-Df_j[1:Gap,]
    N_obtained<-End
  }
  colnames(Matrix_full)<-CNAMES
  return(Matrix_full)
}
TexmexMCarlo_Extrapol<-function (nSample, mexList, pqu_extrapol,mult = 10) 
{
  d <- length(mexList)
  data <- mexList[[1]]$margins$data
  margins <- mexList[[1]]$dependence$margins
  nData <- dim(data)[1]
  which <- sample(1:nData, size = nSample, replace = TRUE)
  MCsampleOriginal <- data[which, ]
  dataLaplace <- Craft_MEX_transform(mexList[[1]]$margins, margins = margins, 
                              method = "mixture")$transformed
  MCsampleLaplace <- dataLaplace[which, ]
  ### Where does the maximum appears ? 
  whichMax <- apply(MCsampleLaplace, 1, which.max) 
  
  ### with the traditional package function
  ### dth <- sapply(mexList, function(l) l$dependence$dth)
  ### with the traditional package function
  ### dqu <- sapply(mexList, function(l) l$dependence$dqu)
  
  ### As we do not use the CEV models with the modeling thresholds,
  ### we need to provide these thresholds 
  dth_extrapol<-sapply(1:d,function(index_coord){
    Lev_quant<-pqu_extrapol[index_coord]
    return(as.numeric(quantile(dataLaplace[, index_coord], 
                               Lev_quant)))
  })
  
  ### When does the maximum is above the desired threshold--> CEV can be used
  whichMaxAboveThresh <- sapply(1:nSample, function(i) MCsampleLaplace[i, 
                                    whichMax[i]] >= dth_extrapol[whichMax[i]])
  mexKeep <- lapply(1:d, function(i) {
    ### to extrapolate--> use new argument pqu_extrapol
    mc <- predict(mexList[[i]], pqu = pqu_extrapol[i], nsim = nSample * 
                        d * mult,smoothZdistribution=TRUE)
    ### keep in memory only the simulations for which the maximum appears 
    ### for the ith model
    mc$data$simulated[mc$data$CondLargest, order(c(i, c(1:d)[-i]))]
  })
  N_computed<-sum(sapply(mexKeep,FUN = function(x){
    return(nrow(x))
  }))
  N_totake<-min(N_computed,nSample)
  nR <- rep(0, d)
  names(nR) <- names(data)
  list_reals<-list()
  cond_max<-whichMaxAboveThresh
  for (i in 1:d) {
    ### When does the simulations from the ith model can be used --> 1) when the ith coordinate is the maximum,
    ### and 2) this maximum is above the threshold => what we are looking for
    replace <- whichMax == i & whichMaxAboveThresh
    cond_i<-replace
    ### instead of using the indexes of replace (much too few),
    ### we estimate the mixing weight--> draw the convenient number for each 
    ### model
    prob_cli<-sum(cond_i)/sum(cond_max)
    print(prob_cli)
    # proportion of simulated samples--> number is N_totake,
    # not automatically nSample
    nReplace_i<-round(N_totake*prob_cli)
    print("To take")
    print(N_totake)
    if (nReplace_i > 0) {
      nR[i] <- nReplace_i
      Mat_cl_i<-as.matrix(mexKeep[[i]])
      Number_sampled_i<-min(nrow(Mat_cl_i),nReplace_i)
      realisations_i<-Mat_cl_i[c(1:Number_sampled_i), ]
      ### Originally the data resampled
      MCsampleOriginal[c(1:Number_sampled_i), ] <-realisations_i
      print(paste0("i=",i))
      print(Number_sampled_i)
      ### Must be after the conditional simulations in our case
      ### If we cannot obtain many samples--> we have a large number of observations
      list_reals[[i]]<-realisations_i
    }
  }
  res <- list(nR = nR, MCsample = MCsampleOriginal, whichMax = whichMax, 
              whichMaxAboveThresh = whichMaxAboveThresh)
  oldClass(res) <- "mexMC"
  return(list("orig_reyy"=res,"list_craft"=list_reals))
}

Craft_MEX_transform<-function (x, margins, r = NULL, method = "mixture", divisor = "n+1", 
                               na.rm = TRUE) 
{
  if (!is.element(method, c("mixture", "empirical"))) 
    stop("method should be either 'mixture' or 'empirical'")
  if (!is.element(divisor, c("n", "n+1"))) 
    stop("divisor can be 'n' or 'n+1'")
  if (is.null(r)) {
    r <- x
    r$transData <- lapply(1:dim(x$data)[2], function(i) x$data[, 
                                                               i])
  }
  transFun <- function(i, x, r, mod, th, divisor, method) {
    x <- x[, i]
    r <- r[[i]]
    mod <- mod[[i]]
    th <- th[i]
    if (divisor == "n") 
      divisor <- length(r)
    else if (divisor == "n+1") 
      divisor <- length(r) + 1
    ox <- order(x)
    r <- sort(r)
    run <- rle(r)
    p <- cumsum(run$lengths)/divisor
    p <- rep(p, run$lengths)
    Femp <- p[sapply(x, function(y) which.min(abs(r - y)))]
    if (method == "mixture") {
      sigma <- exp(mod$coefficients[1])
      xi <- mod$coefficients[2]
      Para <- (1 + xi * (x - th)/sigma)^(-1/xi)
      Para <- 1 - mean(r > th) * Para
      res <- ifelse(x <= th, Femp, Para)
    }
    else res <- Femp
    res[ox] <- sort(res)
    res
  }
  res <- sapply(1:ncol(x$data), transFun, x = x$data, r = r$transData, 
                mod = r$models, th = r$mth, divisor = divisor, method = method)
  dimnames(res) <- list(NULL, names(r$models))
  x$transformed <- margins$p2q(res)
  invisible(x)
}

HTAWN_censorshisp<-function(model_z,prop_z,thresh_censor,ind_ref,z,
                            l_name){
    accept<-FALSE
    while(accept==FALSE){
      Sims_Cond_z<-predict(object= model_z,pqu=prop_z,
                     nsim =2,which = z,
                     smoothZdistribution=TRUE)$data$simulated
      ###Select the first one (breaks if nsim=1)
      Sims_Cond_z<-unlist(Sims_Cond_z[1,])
      d<-length(Sims_Cond_z)
      if(z==ind_ref){
        #c(ind_ref
        return(Sims_Cond_z)
      }else{
        z_minus<-z-1
        Sims_Cond_z<-Sims_Cond_z[l_name]
        Rest_vector<-Sims_Cond_z[1:z_minus]
        M<-max(Rest_vector)
        ### If we use the third model, 
        ### we must not simulate 
        ### extreme obs for the first second obs
        if(M<thresh_censor){
          #c(z,
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
F_find_Quantile<-function(Vect,Target){
  Dce_Target_Tau<-function(Tau){
    Q_<-1-Tau
    Threshold_estim<-as.numeric(quantile(Vect,Q_))
    return(Threshold_estim-Target)
  }
  Qsurv_star_found<-uniroot(Dce_Target_Tau, lower = 0, upper = 1, tol = .Machine$double.eps^0.5)$root
  return(Qsurv_star_found)
}
