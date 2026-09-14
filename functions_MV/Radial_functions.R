### Convergence of the defined "angles" of the radial vector ------------
#################

#' Title
#'
#' @param Df_Rad_Orig 
#' @param l_name 
#' @param k 
#' @param vect_lg 
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param Df_Rad_Orig 
#' @param l_name 
#' @param vector_k 
#' @param Ref_RiskF 
#' @param CPU_hearts 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param Mu_vector 
#' @param Cov_mat 
#' @param g 
#' @param Th_g 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param Mu_vector 
#' @param Cov_mat 
#' @param g 
#' @param Th_g 
#'
#' @return
#' @export
#'
#' @examples
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
### (CEV model) function for threshold selection ---------------

#' Title
#'
#' @param mqu 
#' @param dqu 
#' @param vect_l 
#' @param ind_ref 
#'
#' @return
#' @export
#'
#' @examples
Params_HTawn_one_dqu<-function(mqu,dqu,vect_l,ind_ref){
  
  model_texmex<-texmex::mex(vect_l,mqu = mqu,
                            dqu = dqu,which = ind_ref)
  Dep<-model_texmex$dependence
  Theta<-Dep$coefficients
  Z<-Dep$Z
  n<-nrow(Z)
  d<-ncol(Theta)
  Ref_var<- seq(dqu, 1 - 1/n, length = n)
  Result_per_column<-NA
  Result_per_column2<-NA
  if(d>1){
    ### Verify if test doable
    Sum_ind<-sum(as.numeric(!is.na(Theta[1,])))
    n<-nrow(Z)
    if(Sum_ind==d){
      Resid<-Z
      Absolute_values<-abs(Resid-colMeans(Resid))
      Result_per_column<-apply(Resid,MARGIN = 2,
             FUN = function(x){
               result<-round(cor.test(x,Ref_var,method="kendall")$p.value,
                             3)
               return(result)
             })
      Result_per_column2<-apply(X = Absolute_values,MARGIN = 2,
            FUN = function(x){
              result<-round(cor.test(x,Ref_var,method="kendall")$p.value,
                            3)
              return(result)
            })
      
    }
  }else{
    ### Verify if test doable
    Answer<-!is.na(as.numeric(Theta[1,]))
    if(Answer==TRUE){
      Resid<-as.numeric(Z)
      Absolute_values<-abs(Resid-mean(Resid))
      Result_per_column<-round(cor.test(Resid,
                                        Ref_var, method="kendall")$p.value,3)
      Result_per_column2<-round(cor.test(Absolute_values,
                                         Ref_var,method="kendall")$p.value,3)
      
      
    }
  }
  return(list("Test_1"=Result_per_column,
              "Test_2"=Result_per_column2,
              "Theta"=Theta))
}
#' Title
#'
#' @param Theta_opt 
#' @param vect_name_variables 
#'
#' @return
#' @export
#'
#' @examples
Analysis_diag_HTawn_Conv_Hull<-function(Theta_opt,vect_name_variables){
  
  d<-length(Theta_opt)
  list_GG_for_hull<-list()
  vect_name_corrected<-sapply(vect_name_variables,FUN = Fct_correct_name,
                              replacement="Surge",target="Surcote")
  for(Z in c(1:d)){
    HTawn_z<-Theta_opt[[Z]]
    Modelboot_z<-texmex::bootmex(x = HTawn_z)

    Results<-lapply(Modelboot_z$boot,function(x){
      y<-x$dependence
      matrix<-as.data.frame(t(y[1:2,]))
      colnames(matrix)<-c("a","b")
      rownames(matrix)<-colnames(y)
      df <- tibble::rownames_to_column(matrix
                                       , "variable")
      df$model<-rep(as.character(Z),
                    nrow(df))
      return(df)
    })
    list_GG_for_hull[[Z]]<-do.call(rbind.data.frame,
                                   Results)
  }
  Combinaison_for_hull<-do.call(rbind.data.frame,
                                list_GG_for_hull)
  Combinaison_for_hull$variable<-sapply(Combinaison_for_hull$variable,
            FUN = function(x){
              return(Fct_correct_name(x = x,target = "Surcote",
                                      replacement = "Surge"))
            })
  Vector_shown_chull<-sapply(vect_name_variables,
           FUN = function(x){
             return(Fct_correct_name(x = x,
                                     target = "Surcote",
                                     replacement = "Surge"))
           })
  DEFAULT_colors<-scales::hue_pal()(
    length(Vector_shown_chull))
  Real_labels <- parse(
    text = paste0("u[", Vector_shown_chull, "]")
  )
  Chull<-Combinaison_for_hull %>% 
    group_by(model,variable)%>%  
    slice(chull(a, b))
  
  GG_hull<-ggplot2::ggplot(Combinaison_for_hull,
                           aes(a, b)) +
    geom_point() +
    facet_wrap(~factor(variable,
                       levels = vect_name_corrected))+
    geom_polygon(data =Chull,
                 alpha = 0.3,aes(fill=model)) +
    xlim(0,1)+ylim(-1,1)+
    theme_minimal()+labs(fill="Legend")+
    theme(axis.title=element_text(size=25),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          axis.text = element_text(size=12),
          strip.text = element_text(size = 13),
          legend.direction = "horizontal",
          legend.position = "bottom")
  GG_hull<-GG_hull+
    scale_fill_manual(values=DEFAULT_colors,
                      labels=Real_labels)
  return(GG_hull)
}



#' Title
#'
#' @param vect_dqu 
#' @param Vect_obs 
#' @param MQU 
#' @param chosen_dqu 
#'
#' @return
#' @export
#'
#' @examples
Analysis_diag_HTawn_evol_DQU<-function(vect_dqu,
                                       Vect_obs,MQU,
                                       chosen_dqu){
  
  LNAME<-colnames(Vect_obs)
  LNAME_corrected<-sapply(LNAME,
                          FUN =function(x){
                            Fct_correct_name(x = x,target = "Surcote",
                                             replacement = "Surge")})
  d<-ncol(Vect_obs)
  GG_theta<-ggplot2::ggplot()
  Rult_simplified<-c()
  Indep_whole<-c()
  for(cond_var in c(1:d)){
    Results_Indep_Theta<-sapply(vect_dqu,
                                Params_HTawn_one_dqu,
                                mqu=MQU,vect_l=Vect_obs,
                                ind_ref=cond_var)
    
    Inds_for_analyse<-apply(X = Results_Indep_Theta,
                            FUN = function(x){
                              Ind_nafalse<-sum(rowSums(!is.na(x$Theta)))
                              dim<-ncol(x$Theta)*nrow(x$Theta)
                              return(Ind_nafalse==dim)
                            },MARGIN = 2)
    ThetaI_simplified<-Results_Indep_Theta[,Inds_for_analyse]
    Sub_DQU<-vect_dqu[Inds_for_analyse]
    Result_Independence_test<-t(ThetaI_simplified[c(1:2),])
    Result_several_theta<-do.call(rbind,ThetaI_simplified[3,])
    
    ### Drop NA results
    
    Rult<-melt(Result_several_theta)
    Sub_vars<-c("a","b")
    Nb_rep_DQU<-length(Sub_vars)
    Inds_chosen<-which(Rult$Var1%in%Sub_vars)
    Rult2<-Rult[Inds_chosen,]
    N_rep<-nrow(Rult2)/(Nb_rep_DQU*length(Sub_DQU))
    colnames(Rult2)<-c("param","variable","value")
    chosen_dquj<-chosen_dqu[cond_var]
    Rult2$dqu<-rep(chosen_dquj,
                   nrow(Rult2))
    Rult2$model<-rep(cond_var,
                     nrow(Rult2))
    Rult2$quantile<-rep(c(sapply(Sub_DQU,
                                 FUN = function(x){rep(x,Nb_rep_DQU)})),
                        N_rep)
    Rult_simplified<-rbind(Rult_simplified,
                           Rult2)
    
    ### Independence test
    ######################
    Melting_indep<-melt(apply(Result_Independence_test,
                              MARGIN = 2,FUN = unlist))
    colnames(Melting_indep)<-c("variable","test_used",
                               "pval")
    Melting_indep$test_used<-ifelse(Melting_indep$test_used=="Test_1",
                                    yes = "Z vs Y",
                                    no = "|Z-mean(Z)| vs Y")
    Melting_indep$quantile<-rep(c(sapply(Sub_DQU,
                                         FUN = function(x){rep(x,N_rep)})),
                                Nb_rep_DQU)
    Melting_indep$variable_cond<-rep(cond_var,
                                     nrow(Melting_indep))
    Melting_indep$dqu<-rep(chosen_dquj,
                           nrow(Melting_indep))
    Indep_whole<-rbind(Indep_whole,
                       Melting_indep)
  }
  Dmodel<-length(unique(Rult_simplified$model))
  Real_test<-expression(u[j])
  if(Dmodel==1){
    Real_label<-as.expression(do.call("substitute", 
                                      list(Real_test[[1]], 
                                           list(j=1))))
    Rult_simplified$model<-as.character(Rult_simplified$model)
    Rult_simplified$variable<-sapply(Rult_simplified$variable,
                                     FUN =function(x){
                                       Fct_correct_name(x = x,target = "Surcote",
                                                        replacement = "Surge")})
    GG_theta<-ggplot2::ggplot(Rult_simplified
                              ,aes(x=quantile,y=value,
                                   group=interaction(variable,param,model),
                                   col=param))+
      geom_line()+
      ylab("Estimator")+
      xlab("Dqu")+facet_wrap(~variable,scales = "free_y")+
      labs(col="Legend",
           linetype="Model")+
      geom_vline(aes(xintercept=dqu,
                     linetype="dep_quantile"))+
      scale_linetype_manual(values=2,
                            labels=Real_label)
    
    Indep_whole$variable_cond<-as.character(Indep_whole$variable_cond)
    
    Indep_whole$variable<-sapply(Indep_whole$variable,
                                 FUN =function(x){
                                   Fct_correct_name(x = x,target = "Surcote",
                                                    replacement = "Surge")})
    GG_Indep<-ggplot2::ggplot(Indep_whole,
                              aes(x=quantile,y=pval,
                                  group=interaction(variable_cond,
                                                    test_used,variable),
                                  col=test_used))+
      geom_line()+facet_wrap(~variable)+ylab("p value")+
      xlab("Dqu")+ labs(col="Legend",linetype="Model")+
      geom_hline(yintercept = 0.05,col="black")+
      geom_vline(aes(xintercept = dqu,
                     linetype="dep_quantile"))+
      scale_linetype_manual(values=2,
                            labels=Real_label)
  }else{
    GG_theta<-list()
    GG_Indep<-GG_theta
    for(J in c(1:Dmodel)){
      Real_label<-as.expression(do.call("substitute", 
                                        list(Real_test[[1]], 
                                             list(j=LNAME_corrected[J]))))
      Rult_simplified_i<-Rult_simplified[which(Rult_simplified$model==J),]
      Rult_simplified_i$model<-as.character(Rult_simplified_i$model)
      Rult_simplified_i$variable<-as.character(Rult_simplified_i$variable)
      Rult_simplified_i$variable<-sapply(Rult_simplified_i$variable,
                                         FUN =function(x){
                                           Fct_correct_name(x = x,target = "Surcote",
                                                            replacement = "Surge")})
      
      GG_theta_J<-ggplot2::ggplot(Rult_simplified_i
                                  ,aes(x=quantile,y=value,
                                       group=interaction(variable,param,model),
                                       col=param))+
        geom_line()+
        ylab("Estimator")+
        xlab("Dqu")+
        facet_wrap(~variable,
                   scales = "free_y")+
        #geom_point(aes(shape=model))+
        labs(col="Legend","linetype"="Model")+
        geom_vline(aes(xintercept=dqu,
                       linetype="dep_quantile"))+
        scale_linetype_manual(values=2,
                              labels=Real_label)+
        theme(axis.title = element_text(size=15),
              strip.text=element_text(size=13),
              legend.title = element_text(size=14),
              legend.text=element_text(size=13),
              legend.direction = "horizontal",
              legend.position = "bottom")+
        ggguides::legend_order_guides(linetype = 1, colour = 2)
      GG_theta[[J]]<-GG_theta_J
      
      Indep_whole_i<-Indep_whole[which(Indep_whole$variable_cond==J),]
      Indep_whole_i$variable_cond<-as.character(Indep_whole_i$variable_cond)
      Indep_whole_i$variable<-as.character(Indep_whole_i$variable)
      Indep_whole_i$variable<-sapply(Indep_whole_i$variable,
                                     FUN =function(x){
                                       Fct_correct_name(x = x,target = "Surcote",
                                                        replacement = "Surge")})
      GG_Indep_J<-ggplot2::ggplot(Indep_whole_i,
                                  aes(x=quantile,y=pval,
                                      group=interaction(variable_cond,
                                                        test_used,variable),
                                      col=test_used))+
        geom_line()+
        facet_wrap(~variable)+
        #geom_point(aes(shape=variable_cond))+
        ylab("p value")+
        xlab("Dqu")+
        labs(col="Legend","linetype"="Model")+
        geom_hline(yintercept = 0.05,col="black")+
        geom_vline(aes(xintercept = dqu,
                       linetype="dep_quantile"))+
        scale_linetype_manual(values=2,
                              labels=Real_label)+
        theme(legend.direction = "horizontal",
              legend.position = "bottom",
              axis.title = element_text(size=15),
              strip.text=element_text(size=13),
              legend.title = element_text(size=14),
              legend.text=element_text(size=13))+
        ggguides::legend_order_guides(linetype = 1, colour = 2)
      GG_Indep[[J]]<-GG_Indep_J
    }
  }
  
  return(list("theta"=GG_theta,
              "indep"=GG_Indep))
}

#' Title
#'
#' @param model_z 
#' @param prop_z 
#' @param thresh_censor 
#' @param ind_ref 
#' @param z 
#' @param l_name 
#'
#' @return
#' @export
#'
#' @examples
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
