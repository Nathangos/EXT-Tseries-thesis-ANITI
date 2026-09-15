#' Replace in the input x a substring by another one if 
#' it is present.
#'
#' @param x str. Vector to modify.
#' @param target string. Element to
#' replace in the vector.
#' @param replacement string. Element to use
#' as a replacement.
#'
#' @return
#' @export
#'
#' @examples
Fct_correct_name<-function(x,target,replacement){
  

  b<-as.numeric(str_locate(x,target))
  S<-sum((is.na(b)))
  if(S>=1){
    return(x)
  }else{
    Loc<-substr(x,start = b[1],stop = b[2])
    x_new<-str_replace(x,
                       Loc, 
                       replacement)
    return(x_new)
  }
}
#' new implementation of the plot Acf to increase
#' the the character size in the base plot Acf
#'
#' @param x 
#' @param ci 
#' @param type 
#' @param xlab 
#' @param ylab 
#' @param ylim 
#' @param main 
#' @param ci.col 
#' @param ci.type 
#' @param max.mfrow 
#' @param ask 
#' @param mar 
#' @param oma 
#' @param mgp 
#' @param xpd 
#' @param cex.main 
#' @param verbose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_Acf_modif<-function (x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL, 
                          ylim = NULL, main = NULL, ci.col = "blue", ci.type = c("white", 
                                                                                 "ma"), max.mfrow = 6, ask = Npgs > 1 && dev.interactive(), 
                          mar = if (nser > 2) c(3, 2, 2, 0.8) else par("mar"), oma = if (nser > 
                                                                                         2) c(1, 1.2, 1, 1) else par("oma"), mgp = if (nser > 
                                                                                                                                       2) c(1.5, 0.6, 0) else par("mgp"), xpd = par("xpd"), 
                          cex.main = if (nser > 2) 1 else par("cex.main"), verbose = getOption("verbose"), 
                          ...){
  
  ### plot.acf with the option to add CCF and ACF at diagonal 
  ### and CCF for the anti diagonal
  ci.type <- match.arg(ci.type)
  if ((nser <- ncol(x$lag)) < 1L) 
    stop("x$lag must have at least 1 column")
  if (is.null(ylab)) 
    ylab <- switch(x$type, correlation = "ACF", covariance = "ACF (cov)", 
                   partial = "Partial ACF")
  snames <- x$snames
  with.ci <- ci > 0 && x$type != "covariance"
  with.ci.ma <- with.ci && ci.type == "ma" && x$type == "correlation"
  if (with.ci.ma && x$lag[1L, 1L, 1L] != 0L) {
    warning("can use ci.type=\"ma\" only if first lag is 0")
    with.ci.ma <- FALSE
  }
  clim0 <- if (with.ci) 
    qnorm((1 + ci)/2)/sqrt(x$n.used)
  else c(0, 0)
  Npgs <- 1L
  nr <- nser
  if (nser > 1L) {
    sn.abbr <- x$snames
    if (nser > max.mfrow) {
      Npgs <- ceiling(nser/max.mfrow)
      nr <- ceiling(nser/Npgs)
    }
    opar <- par(mfrow = rep(nr, 2L), mar = mar, oma = oma, 
                mgp = mgp, ask = ask, xpd = xpd, cex.main = cex.main)
    on.exit(par(opar))
    if (verbose) {
      message("par(*) : ", appendLF = FALSE, domain = NA)
      str(par("mfrow", "cex", "cex.main", "cex.axis", "cex.lab", 
              "cex.sub"))
    }
  }
  if (is.null(ylim)) {
    ylim <- range(x$acf[, 1L:nser, 1L:nser], na.rm = TRUE)
    if (with.ci) 
      ylim <- range(c(-clim0, clim0, ylim))
    if (with.ci.ma) {
      for (i in 1L:nser) {
        clim <- clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1, 
                                                   i, i]^2)))
        ylim <- range(c(-clim, clim, ylim))
      }
    }
  }
  for (I in 1L:Npgs) for (J in 1L:Npgs) {
    dev.hold()
    iind <- (I - 1) * nr + 1L:nr
    jind <- (J - 1) * nr + 1L:nr
    if (verbose) 
      message(gettextf("Page [%d,%d]: i =%s; j =%s", I, 
                       J, paste(iind, collapse = ","), paste(jind, collapse = ",")), 
              domain = NA)
    for (i in iind) for (j in jind) if (max(i, j) > nser) {
      frame()
      box(col = "light gray")
    }
    else {
      if(i==j){
        ylab<-"ACF"
      }else{
        ylab<-"CCF"
      } 
      if(x$type=="partial"){
        ylab<-""
      }
      
      clim <- if (with.ci.ma && i == j) 
        clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1, i, j]^2)))
      else clim0
      plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab, 
           ylab = ylab, ylim = ylim, ...)
      abline(h = 0)
      if (with.ci && ci.type == "white") 
        abline(h = c(clim, -clim), col = ci.col, lty = 2)
      else if (with.ci.ma && i == j) {
        clim <- clim[-length(clim)]
        lines(x$lag[-1, i, j], clim, col = ci.col, lty = 2)
        lines(x$lag[-1, i, j], -clim, col = ci.col, lty = 2)
      }
      Name_title<-ifelse(i == j,snames[i],
                         paste(sn.abbr[i], "&", sn.abbr[j]))
      title(Name_title, line = if (nser > 
                                   2) 
        1
        else 2)
    }
    if (Npgs > 1) {
      mtext(paste("[", I, ",", J, "]"), side = 1, line = -0.2, 
            adj = 1, col = "dark gray", cex = 1, outer = TRUE)
    }
    
    dev.flush()
    
  }
  invisible()
}
Create_list_fromALL<-function(l_name,obj_source){
  list_toreturn<-list()
  for(j in c(1:length(l_name))){
    if(typeof(obj_source)=="list"){
      list_toreturn[[l_name[j]]]<-obj_source[[j]]
    }
    if((typeof(obj_source)=="double")&&
       (!is.null(dim(obj_source)))){
      list_toreturn[[l_name[j]]]<-obj_source[,j]
    }
    if((typeof(obj_source)=="double")&&
       (is.null(dim(obj_source)))){
      list_toreturn[[l_name[j]]]<-obj_source[j]
    }
    
  }
  return(list_toreturn)
} 
#' choose the threshold minimising the EQD metric 
#' at a given time.
#'
#' @param var_t vector[float]. Univariate
#' values at each time step.
#' @param nb_threshs int. Number of candidate threshold.
#' @param q_min float. Minimum quantile level.  
#' @param q_max float. Maximum quantile level.  
#' @param plot_graph boolean (FALSE by default). 
#' Ask the printing of the criterion plot
#'
#' @return list. Selected threshold and 
#' ordered quantiles using the criterion.
#' @export
#'
#' @examples
Choice_automatic_thresh_per_variable_time<-function(var_t,nb_threshs,
                                                    q_min,
                                                    q_max,
                                                    plot_graph){
  
  #Original paper
  #qmin=0
  #qmax=0.95
  #nb_threshd=20
  Qtiles_candidates<-seq.int(q_min,q_max,
                   length.out=nb_threshs)
  thresholds1 <- quantile(var_t,Qtiles_candidates)
  Metric_result<-eqd(var_t, 
                     thresh = thresholds1)
  if(plot_graph){
    plot(thresholds1, Metric_result$dists, 
         xlab="Threshold", ylab="Metric value")
  }
  index_used<-which(Metric_result$thresh==thresholds1)
  ### Order thresholds according to the distance
  ### from the best to the worst.
  Order_metric_for_th<-order(Metric_result$dists,
                             decreasing = FALSE)
  Ordered_qtiles_metric_derived<-Qtiles_candidates[Order_metric_for_th]
  return(list("p_u_opt"=Qtiles_candidates[index_used],
              "ordered_quantiles"=Ordered_qtiles_metric_derived)
         )
}

#' Choose for every time step t the criterion
#' minimising the EQD metric.
#' 
#' @param df dataframe[float]. Univariate time series
#' @param nb_threshs int. Number of candidate threshold.
#' @param Qmin float. Minimum quantile level.
#' @param Qmax float. Maximum quantile level.  
#' @param hearts Object from parallel package 
#' to use parallel computations
#'
#' @return
#' @export
#'
#' @examples
Choice_automatic_thresh_per_variable<-function(df,
                      nb_threshs,Qmin, Qmax,hearts){
  
  return(parallel::parApply(cl = hearts,
          X = df,MARGIN = 2,FUN =
          Choice_automatic_thresh_per_variable_time,
          nb_threshs=nb_threshs,q_min = Qmin,
          q_max = Qmax,plot_graph=TRUE))
}

#' Convert every data margin to a heavy-tailed scale
#' (Frechet or GPD) by estimatin the cdf with a
#' mixture approach (empirical for the bulk and GPD for the tail)
#'
#' @param CPU_hearts object from the parallel package. 
#' Used to apply the parallel computation.
#' @param lien_racine string. Link for 
#' the import for the time series measures.
#' @param list_names vector[str]. Vector of variable names.
#' @param file_dates string. Link for the import
#' of date vector.
#' @param n.dens int. Parameter of the 
#' KdE estimate of the bulk of the marginal distributions
#' @param opt_Frech Boolean. Type of marginal 
#' transformation (Frechet or GPD)
#'
#' @return
#' @export
#'
#' @examples
MarTransfo_TS_exts_Mixture<-function(CPU_hearts,lien_racine,list_names,
                                     file_dates,n.dens,
                                     opt_Frech,Nb_Threshs){
  p_U<-list()
  l_Orig<-list()
  j<-1
  l_AD<-list()
  l_Kdens<-list()
  l_ev_theta<-list()
  d<-length(list_names)
  l_Transf<-list()
  for(name_variable in list_names){
    ### Add d to prevent contamination problems.
    lien_donnees<-paste0(lien_racine,name_variable,"_residuals_",
                         d,".csv")
    Donnes<-read.csv(file=lien_donnees)
    l_Orig[[name_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    ### Choice of threshold using metric
    #######
    ### Here, we focus on positive exceedances
    SUM_pos<-apply(Donnes>0,MARGIN = 1,FUN = sum)
    IND_posj<-which(SUM_pos==ncol(Donnes))
    Donnes_pos<-Donnes[IND_posj,]
    Proba_pos<-apply(X = Donnes,MARGIN = 2,
                     FUN = function(x){
                       return(mean(as.numeric(x>0)))
                     })
    p_U_variable<-Choice_automatic_thresh_per_variable(df = Donnes_pos,
                    nb_threshs = Nb_Threshs,Qmin=0, 
                    Qmax = 0.95,hearts=CPU_hearts)
    ### Obtain actual
    p_uOpt<-sapply(c(1:ncol(Donnes)),
              function(k){
                pu_k<-p_U_variable[[k]][["p_u_opt"]]
                proba_posk<-Proba_pos[k]
                proba_pop_pos<-proba_posk*pu_k
                proba_pop_other<-(1-proba_posk)
                cplmt_pu<-proba_pop_pos+proba_pop_other
                ### p_U= weight put on the tail in Opitz.
                return(1-cplmt_pu)
                })

    p_U[[name_variable]]<-p_uOpt
    Vecteur_DIM<-c(1:ncol(Donnes))
    ### Own code
    Resultat_P<-as.data.frame(t(sapply(Vecteur_DIM,
                  f_marginales_all_Pareto,data_to_tf=Donnes,
                     p_u=p_uOpt,
                    n.dens=n.dens)))

    
    PARETO<-Resultat_P$obs
    ### Unit Pareto
    Vtransf<-do.call(cbind.data.frame,PARETO)
    Vunif<-(1-(1/Vtransf))
    if(opt_Frech==TRUE){
      ### Frechet(1)
      Vtransf<--1/log(Vunif)
    }else{
      ### GPD(sig=1,gam=1,u=0)
      Vtransf<-(1/(1-Vunif))-1
    }
   
    colnames(Vtransf)<-substr(colnames(Donnes),
                              start = 2,stop = 3)
    Times_found<-colnames(Vtransf)
    Intermed_list<-list("gamma"=Resultat_P$gamma,
                       "scale"=Resultat_P$scale,
                       "threshold"=Resultat_P$threshold,
                       "p_u"=Resultat_P$p_u,
                       "time_analysed"=Times_found)
    Df_theta<-data.frame(lapply(Intermed_list,
                                        FUN = unlist))
    rownames(Df_theta)<-c(1:nrow(Df_theta))
    l_ev_theta[[name_variable]]<-Df_theta
    
    ### Kdens export
    KDENS<-Resultat_P$kernel_dens
    L<-length(Times_found)
    LK<-lapply(c(1:L),function(t){
      var_x<-KDENS[[t]]
      Intermed_df<-data.frame(var_x$x,var_x$y)
      colnames(Intermed_df)<-c("x","y")
      Nt<-nrow(Intermed_df)
      ### Keep in memory the corresponding time step.
      real_time<-Times_found[t]
      Intermed_df$time_analysed<-rep(real_time,Nt)
      return(Intermed_df)
    })
    l_Kdens[[name_variable]]<-do.call(rbind.data.frame,
                                      LK)
    
    if(name_variable==list_names[1]){
      df<-matrix(NA,ncol=length(list_names),
                 nrow=nrow(Vtransf))
    }
    df[,j]<-apply(X = Vtransf,MARGIN = 1,FUN = calcul_norm_L2)
    j<-j+1
    l_Transf[[name_variable]]<-Vtransf
  }
  # Pareto ------------------------------------------------------------------
  return(list("df"=df,
          "evt_theta"=l_ev_theta,
          "k_dens"=l_Kdens,
          "Transf"=l_Transf))
}
#' #' Obtain the return curve 
#' #'
#' #' @param data_obs dataframe[float]. Univariate 
#' #' time step.
#' #' @param vect_pu vector[float]. Proportion of 
#' #' extreme values for each time step. 
#' #'
#' #' @return Return Curve object.
#' #' @export
#' #'
#' #' @examples
#' RC_emp_transf<-function(data_obs,vect_pu){
#'   
#'   Result_conversion_EXP<-lapply(c(1:ncol(data_obs)),function(x){
#'     return(RC_emp_transf_per_time(col_obs = data_obs[[x]],
#'                                   p_u = vect_pu[x]))
#'   })
#'   return(Result_conversion_EXP)
#' }
#' #' RC_emp_transf_per_time
#' #'
#' #' @param col_obs vector[float]. Univariate values at a given time 
#' #' step
#' #' @param p_u float. Proportion of extreme values. 
#' #'
#' #' @return Return Curve object.
#' #' @export
#' #'
#' #' @examples
#' RC_emp_transf_per_time<-function(col_obs,p_u){
#'   
#'   return(ReturnCurves::margtransf(data = col_obs,
#'                                   qmarg = p_u))
#' }

fct_extract_Transf<-function(time_list){
  return(time_list[["unif_convert_t"]])
}
fct_extract_ParamsEGPD<-function(time_list){
  return(time_list[["params_egpd"]])
}
fct_extract_InitEGPD<-function(time_list){
  return(time_list[["Init"]])
}
fct_extract_FittingEGPD<-function(time_list){
  return(time_list[["Model_t"]])
}
#' Represent with a ggplot a sample of multivariate time series for one 
#' or two groups (option_facet_wrap)
#'
#' @param L_nameV 
#' @param Max_time 
#' @param Melting_df 
#' @param option_facet_wrap 
#'
#' @return Ggplot with df from melting operation.
#' @export
#'
#' @examples
Fct_cplmt_ggplot<-function(Times_available,vect_name_variable,Melting_df,
                           option_facet_wrap){
  
  Comment<-"Time(hour) with respect to tidal peak"
  d<-length(vect_name_variable)
  Time_index_MV<-c()
  Time_index_first<-c()
  Stotal<-c()
  BREAKS<-c(1)
  LABELS<-c()
  S_sub<-0
  for(jn in c(1:d)){
    nameV<-vect_name_variable[jn]
    new_TA<-as.numeric(Times_available[[nameV]])
    Number_times_seen<-length(new_TA)
    ### Times of concatenated Omega
    Cplmt<-c(1:Number_times_seen)+S_sub
    Time_index_first<-c(Time_index_first,
                        Cplmt)
    S_sub<-Number_times_seen+S_sub
    
    ### Concatenated real times. 
    L<-length(Time_index_first)
    Stotal<-c(Stotal,Time_index_first[L])
    Time_index<-(new_TA-19)/6
    Time_index_MV<-c(Time_index_MV,
                     Time_index)
    Midway<-which(Time_index==0)
    BREAKS<-c(BREAKS,
              c(Cplmt[Midway],
                Cplmt[length(Cplmt)]))
    ### Labels change according to the position 
    ### of the variable
    t_modif<-Time_index[length(Time_index)]
    if(t_modif<0){t_modif<-ceiling(t_modif)
    }else{t_modif<-floor(t_modif)}
    opt_beg<-c()
    if(jn==1){
      BEGTIME<-as.character(Time_index[1])
      opt_beg<-BEGTIME
    }
    if(jn!=d){
      nameVfollow<-vect_name_variable[jn+1]
      new_TA_follow<-as.numeric(Times_available[[nameVfollow]])
      L_follow<-length(new_TA_follow)
      Scd_part<-(new_TA_follow[1]-19)/6
      t_modif2<-Scd_part
      if(t_modif2<0){t_modif2<-ceiling(t_modif2)
      }else{t_modif2<-floor(t_modif2)}
      ENDTIME<-paste(t_modif,
                     t_modif2)
      Cplmt_label<-c(opt_beg,"0",
                     ENDTIME)
    }
    if(jn==d){
      ENDTIME<-paste0(t_modif)
      Cplmt_label<-c("0",
                     ENDTIME)
    }
    LABELS<-c(LABELS,
              Cplmt_label)
  }
  ### whole construction
  Ogn_indexes<-c(1:length(LABELS))
  ind_final<-length(Ogn_indexes)-1
  final_indexes<-Ogn_indexes[2:ind_final]
  sub_breaks<-BREAKS[2:ind_final]
  Inds_nonpair<-which(final_indexes%%2!=0)
  sub_whole<-sub_breaks[Inds_nonpair]
  vline_full<-sub_whole
  print(vline_full)
  
  ### dash construction
  Value_y_for_spec<-max(Melting_df$value)
  Times_T<-which(LABELS=="0")
  sub_dash<-BREAKS[Times_T]
  vline_dashed<-sub_dash
  vect_name_variable_Corrected<-sapply(vect_name_variable,
                FUN = function(x){
                  return(Fct_correct_name(x=x,
                        target = "Surcote",
                        replacement = "Surge"))
                })
  Df_for_annotate<-cbind.data.frame(sub_dash,
              rep(Value_y_for_spec,length(sub_dash)),
              vect_name_variable_Corrected)
  colnames(Df_for_annotate)<-c("xspec","yspec","label_spec")
  GG_present_object<-ggplot()+
    geom_text(
      data=Df_for_annotate,
      aes(x=xspec,y=yspec,label=label_spec),
      size=12
    )
 
  if("nameFC"%in%colnames(Melting_df)){
    Melting_df$origin<-paste(Melting_df$nameFC,"-",
                             Melting_df$origin)
    GG_present_object<-GG_present_object+
      geom_line(data = Melting_df,aes(x=time,y=value,
            group=interaction(index,
                              origin),
               col=index))
  }else{
    GG_present_object<-GG_present_object+
      geom_line(data = Melting_df,aes(x=time,y=value,
            group=interaction(index),
            col=index))
  }
  if(option_facet_wrap==TRUE){
    GG_present_object<-GG_present_object+
      facet_wrap(~origin)
  }
  GG_present_object<-GG_present_object+
    geom_line()+
    guides(col="none")+
    guides(col="none")+
    labs(col="Legend")+
    ylab("value (-)")+
    xlab(Comment)+
    geom_vline(xintercept=vline_dashed,
               linetype="dashed")+
    geom_vline(xintercept = vline_full)+
    scale_x_continuous(breaks = BREAKS,
                       labels = LABELS)+
    theme(axis.title=element_text(size=20),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          axis.text = element_text(size=12),
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 15))
  return(GG_present_object)
}
#' Use the 
#'
#' @param z: int. Time index. 
#' @param Data_pos: dataframe. Time series of positive marginals. 
#' @param show_EGPD: Bool. Show or not the details of EGPD fitting. 
#' @param list_params_EGPD list[vector,vector,vector]. EGPD parameter obtained
#' at each time step. 
#'
#' @return
#' @export
#'
#' @examples
Convert_time_z<-function(z,Data_pos,show_EGPD,
                         list_params_EGPD){
  print(paste0("Result for time ",z))
  n.cyc<- list_params_EGPD[["n.cyc"]]
  mu.step<- list_params_EGPD[["mu.step"]]
  sigma.step<- list_params_EGPD[["sigma.step"]]
  nu.step<- list_params_EGPD[["nu.step"]]
  tau.step<- list_params_EGPD[["tau.step"]]
  EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")
  Col_t<-Data_pos[,z]
  Col_pos_t<-Col_t[which(Col_t>0)]
  db<-as.data.frame(Col_pos_t)
  colnames(db)<-c("x")
  con <- gamlss.control(n.cyc = n.cyc,mu.step = mu.step,
                        sigma.step = sigma.step, nu.step =nu.step,
                        tau.step = tau.step,autostep=TRUE,
                        trace = show_EGPD)
  con.i<-glim.control(glm.trace = FALSE)

  ## Nu.start initialisation using moments.
  Th_beg<-quantile(x = Col_pos_t,0.05)
  Lower_tail<-Col_pos_t[which(Col_pos_t<Th_beg)]
  Moment_1<-mean(Lower_tail)
  nu.start<-as.numeric((1-(Moment_1/Th_beg))^(-1)-1)
  
  Opt_thresholds<-Choice_automatic_thresh_per_variable_time(var_t = Col_pos_t,
            q_min = 0,q_max = 0.95,
             nb_threshs = 20,plot_graph = TRUE)
  Candidates_ordered<-Opt_thresholds[["ordered_quantiles"]]
  for(Q_opt in Candidates_ordered){
    Th_opt<-as.numeric(quantile(Col_pos_t,Q_opt))
    FIT_pos<-extRemes::fevd(x = Col_pos_t,threshold = Th_opt,
                            method = "Lmoments",
                            type = "GP")$results
    Shape<-FIT_pos[2]
    ## Use gpd property to get the scale at 0.
    Scale<-FIT_pos[1]-Th_opt*Shape
    INIT<-c(Shape,Scale)
    print(INIT)
    ### EGPD optimisation
    result_EGPD_t <- tryCatch(
      {
        Fitting_time_t <- gamlss(x~1, 
           data=db,family = EGPD1Family(mu.link = "identity"),
           control = con,mu.start=INIT[1],
           sigma.start=INIT[2],nu.start=nu.start,
           i.control=con.i,method=CG())
      },
      error = function(e) {
        # Handle the error
        cat("Error occurred:", conditionMessage(e), "\n")
        return(NA)
      }
    )
    if(length(result_EGPD_t)!=1){
      break
    }
  }

  
  muFit <-fitted(result_EGPD_t,"mu")[1]
  sigmaFit <- predict(result_EGPD_t,what="sigma", 
                      type="response")[[1]]
  nuFit <- predict(result_EGPD_t,what="nu", 
                   type="response")[[1]]
  Params_egpd<-list("mu"=muFit,"sigma"=sigmaFit,"nu"=nuFit)
  vect_unif_conver_t<-pEGPDModel1(q = Col_t,mu =muFit,
                                  sigma = sigmaFit,
                                  nu = nuFit)
  return(list("params_egpd"=Params_egpd,
              "unif_convert_t"=vect_unif_conver_t,
              "Model_t"=result_EGPD_t,
              "Init"=c(INIT,nu.start)))
}

#' Title
#'
#' @param CPU_hearts 
#' @param lien_racine 
#' @param list_names 
#' @param opt_Frech 
#' @param n.dens 
#' @param p_U 
#' @param file_dates 
#' @param show_EGPD 
#' @param list_params_EGPD 
#'
#' @return
#' @export
#'
#' @examples
MarTransfo_TS_EXTGP_MV<-function(CPU_hearts,lien_racine,
                                 list_names,opt_Frech,n.dens,p_U,
                                 file_dates,show_EGPD,
                                 list_params_EGPD){
  
  l_ALL<-list()
  l_Orig_all<-list()
  l_Orig<-list()
  l_AD<-list()
  LIST_SUB_orig<-list()
  LIST_SUB_TRANSFO<-list()
  j<-1
  df<-list()
  d<-length(list_names)
  for(name_variable in list_names){
    
    ### ADD d to prevent contamination error
    lien_donnees<-paste0(lien_racine,name_variable,
              "_residuals_",d,".csv")
    Donnes<-read.csv(file=lien_donnees)
    rownames(Donnes)<-c(1:nrow(Donnes))
    l_Orig_all[[name_variable]]<-Donnes
    dates_import<-read.csv(file=file_dates)[,2]
    dates_import<-dates_import[-1]
    d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
    Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
    NPY<-nrow(Donnes)/Nb_annees
    AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,fonction_seuil = p_U[[name_variable]],
                          n.dens = n.dens,
                              name=name_variable,type_entree=type_entree,
                              dates_prises=dates_import,j_show = 19)
    print("passed !")
    l_AD[[name_variable]]<-AN_GPD
    
    SUM_pos<-apply(Donnes>0,MARGIN = 1,FUN = sum)
    IND_posj<-which(SUM_pos==ncol(Donnes))
    if(j==1){
      IND_select<-IND_posj
    }
    IND_select<-intersect(IND_select,IND_posj)
    j<-j+1
  }
  return(l_Orig_all)
  for(name_variable in list_names){
    Donnes<-l_Orig_all[[name_variable]]
    
    rownames(Donnes)<-c(1:nrow(Donnes))
    Vecteur_DIM<-c(1:ncol(Donnes))
    # 1) Use GAM to analyse short-tail distribs ----------------------------------
    Donnes_pos<-Donnes[IND_select,]
    l_Orig[[name_variable]]<-Donnes_pos
    Unif_EGPD<-matrix(NA,nrow = length(IND_select),
                      ncol=ncol(Donnes))
    print("#########")
    print(name_variable)
    ALL_results<-lapply(colnames(Donnes),FUN =Convert_time_z,
                        Data_pos=Donnes_pos,
                        show_EGPD=show_EGPD,
                        list_params_EGPD=list_params_EGPD)
    Unif_EGPD<-cbind.data.frame(lapply(X = ALL_results,
                            FUN = fct_extract_Transf))
    Vect_egpd<-lapply(X = ALL_results,
                      FUN = fct_extract_ParamsEGPD)
    Vect_Init<-lapply(X = ALL_results,
           FUN = fct_extract_InitEGPD)
    Vect_fitting<-lapply(X = ALL_results,
                         FUN = fct_extract_FittingEGPD)
    
    if(opt_Frech){
      FINAL_transfo<--1/log(Unif_EGPD)
    }else{
      FINAL_transfo<-(1/(1-Unif_EGPD))-1
    }
    FINAL_transfo<-as.data.frame(FINAL_transfo)
    colnames(FINAL_transfo)<-colnames(Donnes)
    ### Remove from FINAL_transfo's colnames
    colnames(FINAL_transfo)<-substr(colnames(FINAL_transfo),
                                    start = 2,
                                  stop = 3)
    l_ALL[[name_variable]]<-list("transf"=FINAL_transfo,
                                "params_transfo"=Vect_egpd,
                                "Fitting"=Vect_fitting,
                                "INIT"=Vect_Init)
    df[[name_variable]]<-apply(X = FINAL_transfo,MARGIN = 1,
                  FUN = calcul_norm_L2)

    }
    return(list("resume"=l_ALL,"df"=as.data.frame(df),"GPD"=l_AD,
            "orig"=l_Orig,"orig_all"=l_Orig_all,
            "SUB_set"=IND_select,"type_transfo"="EXTGPD"
            ))
}


#' Simulate extreme multivariate time series using 
#' some results from regular variations.
#'
#' @param result_transformation : result of previous function.
#' @param l_variables : list[str]. List of variables.
#' @param Q_thresh : float. Proportion of extreme events (to choose wisely). 
#' @param NbScores_Omega : int. Number of scores of the PCA. 
#' @param M : int. Number of simulated time series. 
#' @param rotations_available Boolean. If TRUE, copula 
#' models based on the rotations can be used.  
#' @param Params_risk_Function List. List containing several key arguments such 
#' as the definition of the compound risk function, graphical options and 
#' some functions used for optimisation purpose. 
#' @param root_for_export String. Directory path for the exported elements. 
#' @param list_nb_scores list or int. Number of PC scores kept. 
#' @param One_PCA_base Boolean. If True, the PCA is applied 
#' on the concatenaed angles. 
#' @param f_transf Optional function to apply on the angles
#'  before PCA. Not used.
#' @param f_transf_inv Optional function to go back to
#'  original angles. Not used.
#' @param opt_Frech 
#' @param cols_ggplot vector. Graphic options used in ggplot. 
#' @param Common_theme ggplot2::theme object. Argument to 
#' impose the same graphical options on every ggplot2 graphic. 
#'
#' @return List. Extreme multivariate time series 
#' with several modeling details
#' 
#' @export
#'
#' @examples
Simul_MV_residuals<-function(result_transformation,l_variables,Q_thresh,
                             NbScores_Omega,M,rotations_available,
                             Params_risk_Function,root_for_export,
                             list_nb_scores,One_PCA_base,f_transf,
                             f_transf_inv,opt_Frech,
                             cols_ggplot,Common_theme){
  
  
  Shape_evd<-Params_risk_Function[["Shape_parameter"]]
  # CF Kokozka ---------------------------------------------------------------
  L<-300
  vect_k<-c(20:L)
  Array_simul<-array(data = NA,dim = c(length(l_variables),
                M,37))
  
  list_obs_exts<-list()
  Times_per_FCG<-list()
  list_all<-result_transformation$resume
  for(name_FCG in names(list_all)){
    RAW<-colnames(list_all[[name_FCG]]$transf)
    Times_per_FCG[[name_FCG]]<-substr(RAW,start = 2,
                           stop = 3)  
  }
  Vect_l_function<-result_transformation$df
  
  # RiskF -------------------------------------------------------------------
  ### ----------------------------------------------------------------------
  for(ztilde in c(1:ncol(Vect_l_function))){
    Graphics_estimators_gamma(series = Vect_l_function[,ztilde],
          vect_k = vect_k,
          NB_years = 37,
          Title_graphic = paste0("Shape parameter of l(T("
                                 ,l_name[ztilde]
                                 ,")"))
  }
  # Choice of the individuals -----------------------------------------------
  RiskF_data<-apply(Vect_l_function,MARGIN = 1,
        FUN = List_Params_RF[["function"]])
  Threshold<-quantile(RiskF_data,
                         Q_thresh)
  Indices_exts<-which(RiskF_data>Threshold)
  
  N<-nrow(Vect_l_function)
  L_Transf_vect_l<-list()
  L_Convert<-list()
  
  dir_AIC<-"RisKFunctions/evol_AIC_param"
  list_mixt<-list()
  list_Qfound<-list()
  ### the variable Approach corresponds
  ### to the tool used to model
  ### the distribution of the radial components
  
  Approach<-Params_risk_Function[["RF_Approach"]]
  Vect_l_transf<-NA
  Vect_Pareto_margins<-NA
  if(Approach!="HTawn"){
    ### If necessary, approach 
    ### the marginal distributions of 
    ### the radial vector
    for(nameV in l_variables){
      #log
      VECT_lj<-log(Vect_l_function[,nameV]+1)
      # EGPD way
      Values_gp<-mev::gp.fit(VECT_lj,threshold = 0)$est
      Th_beg<-quantile(x = VECT_lj,0.10)
      Lower_tail<-VECT_lj[which(VECT_lj<Th_beg)]
      Moment_1<-mean(Lower_tail)
      Pstart<-as.numeric((1-(Moment_1/Th_beg))^(-1)-1)
      Initialisation<-c(Pstart,Shape_chosen,0)
      PWM_EXTGP<-mev::fit.extgp(data = VECT_lj,model = 1,
                   init=Initialisation,
                       method="pwm",R=20)
      Theta<-PWM_EXTGP$fit$pwm
      print(c(nameV,Theta))
      Name_file_plot<-paste0(root_for_export,"RisKfunctions/",
                  nameV,"_EGPD_fitting_RiskF.png")
      MLE_EXTGP<-mev::fit.extgp(data = VECT_lj,model = 1,
                    init=Theta,method="mle",
                    R=20)$fit$mle
      u_vec_gom <- quantile(VECT_lj,
                            probs = seq(0.2, 0.9, by = 0.05))
      Choice_threshr_j<-threshr::ithresh(data = VECT_lj, 
                       u_vec = u_vec_gom,
                       n_v = 2, prior = "mdi",
                       h_prior = list(a = 0.6))
      
      Sthresh_j<-summary(Choice_threshr_j)
      plot(Choice_threshr_j)
      Q_foundj<-median(Sthresh_j[,4])
      list_Qfound[[nameV]]<-Q_foundj/100
  
      UNIF_extgp<-sapply(VECT_lj,mev::pextgp,kappa=MLE_EXTGP[["kappa"]],
                              sigma = MLE_EXTGP[["sigma"]],
                              xi=MLE_EXTGP[["xi"]])
      
      Quantile_LEVELS<-c(1:length(VECT_lj))/(length(VECT_lj)+1)
  
      # MLE_EXTGP
      L_Transf_vect_l[[nameV]]<-MLE_EXTGP
      Frechet_j<-(-log(UNIF_extgp))**(-1)
      print(paste0("Frechet test for conversion of l(",nameV,")"))
      print(goftest::ad.test(Frechet_j,null = extRemes::"pevd",
                       shape=1,scale=1,loc=1,type="GEV")$p.value)
      L_Convert[[nameV]]<-Frechet_j
    }
  Vect_l_transf<-as.data.frame(L_Convert)
  
  ### Regular variations hypothesis
  Simple_MRV_HRV_analysis(couple_UV =  Vect_l_transf,
                          root_graphics = paste0(root_for_export,
                            "RisKfunctions/"),
                          Vect_k = c(20:500),
                          l_name = l_variables,q = 0.80)
  Vect_lunif<-exp(-Vect_l_transf^(-1))
  Vect_Pareto_margins<-apply(Vect_lunif,
                    MARGIN=2,FUN = evd::qgpd,loc=0,
                              scale=1, shape=1)
  l_Pareto<-apply(Vect_Pareto_margins,MARGIN = 1,
                  FUN = List_Params_RF[["function"]])

  Seuil_Margins_Pareto<-quantile(l_Pareto,
                                 Q_thresh)
  IndExceedances_Pareto<-which(l_Pareto>Seuil_Margins_Pareto)
  Exceedances_Pareto<-l_Pareto[Indices_exts]

  Theta_obtained<-atan(Vect_Pareto_margins[,2]/Vect_Pareto_margins[,1])
  png(filename= paste0(root_for_export,
                       "/RisKfunctions/ang_dens_R_",
                       l_variables[1],"_",
                       l_variables[2],".png"),
      width=600,height=800)
  plot(density(Theta_obtained),
      main="Angular density of the risk functions",
      cex.lab = 1.5,
      cex.main=1.5)
  abline(v=0,col="red")
  abline(v=pi/2,col="red")
  dev.off()
  ### Distribution of angular comp for exceedances
  Exced_PTO_marg<-Vect_Pareto_margins[Indices_exts,]
  l_exts<-l_Pareto[Indices_exts]
  Angle_exceeds<-t(t(Exced_PTO_marg)%*%diag(l_exts^(-1)))

  df_angle_exceeds<-as.data.frame(Angle_exceeds)
  colnames(df_angle_exceeds)<-c("V1","V2")
  GG_dens_Ang<-ggplot(data = df_angle_exceeds,aes(x=V1,y=V2))+
    geom_point(col="blue")+
    geom_xsidedensity(alpha=0.5)+
    geom_ysidedensity(alpha=0.5)+
    xlab(paste0("Angular component for risk of ",
                l_name[1]))+
    ylab(paste0("Angular component for risk of ",
                l_name[2]))
  ggsave(filename = paste0(root_for_export,
         "/RisKfunctions/Angular_val_exts_",
         l_variables[1],"_",
         l_variables[2],".png"),
         plot = GG_dens_Ang,width=6,height=8)
  Theta_opt<-NA

  }
  ## (1) Approximate the distribution of radial tail components --------
  ####################
  if(Approach=="AD"){
    if(Params_risk_Function[["parametric"]]==TRUE){
      Theta_opt<-Estim_param_RF_homogeneous(Params_risk_Function= Params_risk_Function,
                Q_thresh =1-Q_thresh,
                Vect_l_function = Vect_l_transf,
                d = ncol(Vect_l_transf))
      if(Params_risk_Function[["name_RF"]]=="max"){
        Result_AIC<-as.numeric(as.data.frame(Theta_opt$result_AIC))
        Names_model<-names(Theta_opt$result_AIC)
        Result_AIC<-cbind.data.frame(Names_model,
                                     Result_AIC)
        colnames(Result_AIC)<-c("model","AIC_found")
        Result_AIC$index<-c(1:nrow(Result_AIC))
        GG1<-ggplot(data=Result_AIC,aes(x=index,y=AIC_found,col=model))+
          geom_point()+
          labs(col="MEV model")
        ggsave(filename= paste0(root_for_export,
                                dir_AIC,"_",
                                Approach,".png"),
               plot=GG1,width=8,height=6)
        dev.off()
        Theta_opt<-Theta_opt$alpha
      }
    }
  }
  
  if(Approach=="Gauss"){
    Vect_l_norm<-apply(X = Vect_lunif,
                       FUN = function(x){return(qnorm(x))},
                       MARGIN = 2)
    g_Gauss<-apply(Vect_l_norm,MARGIN = 1,
                   FUN = List_Params_RF[["function"]])
    Mu_Gauss<-quantile(g_Gauss,Q_thresh)
    
    Mu_vector<-colMeans(Vect_l_norm)
    Cov_mat<-cov(Vect_l_norm)
    ### plot MV diag plot
    MVN_diag_plot<-MVN::multivariate_diagnostic_plot(Vect_l_norm)
    png(filename = paste0(root_for_export,"RisKfunctions/",
                          l_name[1],"_",l_name[2],"_Maha_distances.png"))
    print(MVN_diag_plot)
    dev.off()
    ### Summary of results
    Asymp_test<-MVN::mvn(Vect_l_norm, mvn_test="mardia",
                   bootstrap = FALSE)
    print(Asymp_test$multivariate_normality)
    print("### MVN test")
    Mardia_test<-round(MVN::mardia(Vect_l_norm)$p.value,2)
    print(Mardia_test)
    Royston_test<-round(MVN::royston(Vect_l_norm)$p.value,2)
    Energy_test<-round(MVN::energy(Vect_l_norm)$p.value,2)
    HZ_test<-round(MVN::hz(Vect_l_norm)$p.value,2)
    Name_models<-c("mardiaSkew","mardiaKurtosis","Royston",
                   "Energy test","HZ")
    Result_models<-c(Mardia_test, Royston_test,
                     Energy_test,HZ_test )
    Table_Gauss_hyp<-cbind.data.frame(Name_models,Result_models)
    colnames(Table_Gauss_hyp)<-c("name_test",
                                 "p_value")
    write.csv(file=paste0(root_for_export,"RiskFunctions/",
                          l_name[1],"_",l_name[2],"_MVN_test_results.csv"),
              x = Table_Gauss_hyp)
  }
  if(Approach=="HTawn"){
    IREF<-1
    Q_from_mixt<-apply(Vect_l_function,
          MARGIN = 2,
          FUN = function(col_){
            estimated_surv_prob<-F_find_Quantile(Vect = col_,
                                Target = Threshold)
            
          })
    DQU_each_scenario<-1-Q_from_mixt
    ### Choose the du used during modelling step.
    # Use the one used for future simulations
    # DQU_modeling_Htawn<-Q_from_mixt
    
    d<-length(l_variables)
    DQU_modeling_Htawn<-rep(NA,d)
    nb_threshs<-50
    Qtile_candidates<-seq.int(0.70,0.98,
              length.out=nb_threshs)

    MQU<-Params_risk_Function[["HTAWN_params_MARG"]]
    DQU_modeling_Htawn<-Params_risk_Function[["HTAWN_params"]]
  
    Graphics_diags<-Analysis_diag_HTawn_evol_DQU(
             vect_dqu = Qtile_candidates,
             chosen_dqu=DQU_modeling_Htawn,
             MQU =MQU,Vect_obs = Vect_l_function)
    
    GG_theta<-Graphics_diags[["theta"]]
    GG_Indep<-Graphics_diags[["indep"]]
    if(is.list(GG_theta)){
      Ensemble<-c(1:length(GG_theta))
      for(elt in Ensemble){
        #Ensemble_rest<-l_name[Ensemble[-elt]]
        GGtheta_elt<-GG_theta[[elt]]+
          Common_theme
        GGIndep_elt<-GG_Indep[[elt]]+
          Common_theme
        ggsave(filename = paste0(root_for_export,"RiskFunctions/d=",
                                 length(l_name),"_evol_HTparams_model",elt,".png"),
               plot = GGtheta_elt,width = 10,height = 6)
        ggsave(filename = paste0(root_for_export,"RiskFunctions/d=",
                                 length(l_name),
                          "_result_HTindep_model",elt,".png"),
               plot =  GGIndep_elt,width = 10,height = 6)
      }
    }else{
      ggsave(filename = paste0(root_for_export,"RiskFunctions/",
                               l_name[1],"_",l_name[2],"_d=",
                               length(l_name),"_evol_HTparams.png"),
             plot = GG_theta,width = 8,height = 6)
      ggsave(filename = paste0(root_for_export,"RiskFunctions/",
                               l_name[1],"_",l_name[2],"_d=",
                               length(l_name),"_result_HTindep.png"),
             plot =  GG_Indep,width = 8,height = 6)
    }
    # HTawn modeling ---------------------------------------------------------
    ############
    Model_Htawn_all<-texmex::mexAll(
            data = Vect_l_function,
            mqu = MQU,
            dqu = DQU_modeling_Htawn)
    Theta_opt<-Model_Htawn_all
    list_GG_cases<-list()
    for(Z in c(1:d)){
      Vars_present<-c(1:d)[-Z]
      Names_other<-sapply(l_variables[Vars_present],
                          FUN = Fct_correct_name,
                          target="Surcote",replacement="Surge")
      D_rest<-length(Vars_present)
      Name_Var_cond<-l_variables[Z]
      L_increment<-1
      ### Graphics to visualize what appears
      for(vars in Vars_present){
        #### Modify ggplot vline variable name
        Var_cond_replaced<-Fct_correct_name(
          x = Name_Var_cond,target = "Surcote",
          replacement = "Surge"
        )
        Other_var<-l_variables[vars]
        ### Var cond first, other var second
        PAIR_for_graphic<-c(Name_Var_cond,Other_var)
        PAIR_for_graphic_corrected<-sapply(PAIR_for_graphic,
                 FUN = function(x){
                   return(Fct_correct_name(x = x,
                                           target = "Surcote",
                                           replacement = "Surge"))
                                           })
        HTawn_z<-Theta_opt[[Z]]
        
        ### Residual against predictor
        model<-HTawn_z$dependence
        MARGINS<-HTawn_z$margins$transformed
        Th<-model$dth
        Var_cond<-MARGINS[,model$conditioningVariable]
        Inds_exts<-which(Var_cond>Th)
        Reg_HT<-Var_cond[Inds_exts]
        Resid<-model$Z
        VARIABLE<-c(sapply(Names_other,
                           FUN=function(x){
                             Category<-x
                             rep(Category,nrow(Resid))}))
        REGRESSOR<-rep(x = Reg_HT,ncol(Resid))
        RESID<-c(Resid)
        DF<-cbind.data.frame(VARIABLE,REGRESSOR,RESID)
        colnames(DF)<-c("variable","reg_","resid_")
        GG_trend_Z_Yreg<-ggplot(DF,aes(x=reg_,y=resid_))+
          facet_wrap(~variable)+
          geom_point()+
          geom_smooth(aes(col="Polynomial regression"),
                      alpha=0.4)+
          labs(col="Legend")+
          xlab(paste0("Conditioning variable (",Var_cond_replaced,")"))+
          ylab("Residual")+
          Common_theme+theme(
                legend.direction = "horizontal",
                legend.position = "bottom")
        ggsave(plot = GG_trend_Z_Yreg,
               filename = paste0(root_for_export,"RiskFunctions/d=",
                  length(l_name),"_trend_Reg_Resid_HT_model",
                  Z,".png"),
               width=10,height=6)

        ### Input from HTawn
        
        Data_tfed<-as.data.frame(MARGINS[,PAIR_for_graphic])
        colnames(Data_tfed)<-sapply(c(1:length(PAIR_for_graphic)),
        function(x){
          return(paste0("V",x))
        })
        Th_used_j<-as.numeric(HTawn_z$dependence$dth)
        Data_tfed$inds_exts_Z<-as.character(
          as.numeric(Data_tfed[,1]>Th_used_j))
        
        Label_1<-latex2exp::TeX(sprintf("$\\ell(T(%s))$",
                                  PAIR_for_graphic_corrected[1]))
        Label_2<-latex2exp::TeX(sprintf("$\\ell(T(%s))$",
                                  PAIR_for_graphic_corrected[2]))
        GG_z_model<-ggplot(data = Data_tfed,aes(x=V1,y=V2,
                                                col=inds_exts_Z))+
          geom_point()
        Real_test<-expression(u[Z])
        Real_label<-as.expression(do.call("substitute", 
              list(Real_test[[1]], 
                   list(Z=Var_cond_replaced))))
        GG_z_model<-GG_z_model+
          geom_vline(aes(xintercept=Th_used_j,
                         linetype="dep_quantile"))+
          scale_linetype_manual(values=2,
                                labels=Real_label)
        GG_z_model<-GG_z_model+ 
          labs(col="Legend",linetype="Model")+
          xlab(Label_1)+
          ylab(Label_2)+
          Common_theme+theme(
                legend.direction = "horizontal",
                legend.position = "bottom")
        if(L_increment==1){
          LGD_input<-ggpubr::get_legend(GG_z_model)
          LGD_HTAWN<-plot_grid(LGD_input)
          GG_z_model<-GG_z_model+
            guides(col="none",linetype="none")
          WHOLE<-GG_z_model
          
        }
        if(!L_increment%in%c(1,D_rest)){
          if(L_increment%%2==0){
            ### 2 graphics per row
            WHOLE<-WHOLE+GG_z_model+
              guides(col="none",linetype="none")
          }else{
            GG_transit<-GG_z_model+
              guides(col="none",linetype="none")
            ### We jump the line if it is not the case
            WHOLE<-WHOLE/GG_transit
          }

        }
        #ADD +plot_layout(guides="collect")
        if(L_increment==D_rest){
          if(L_increment%%2==0){
            GG_transit<-GG_z_model+
              guides(col="none",linetype="none")
              
            WHOLE<-WHOLE+GG_transit
          }else{
            GG_transit<-GG_z_model+
              guides(col="none",linetype="none")
            ### Fill the last row with plot_spacer
            Second_row<-plot_spacer() + GG_transit+ plot_spacer()+
              plot_layout(widths = c(0.5, 1,0.5))
            WHOLE<-WHOLE/Second_row
          }
        }
        ### Increment the index for next graphic to add
        L_increment<-L_increment+1
      }
      WHOLE_RF<-WHOLE+
        plot_layout(axis_titles = "collect_x")
      Base_RF_modif<-plot_grid(WHOLE_RF,
               LGD_HTAWN,nrow=2,rel_heights =c(10,1))
      WHOLE_RF_modif<-ggdraw()+
        draw_plot(Base_RF_modif)+
        theme(plot.background = element_rect(fill = "white", 
                                             color = NA))
      ggsave(plot = WHOLE_RF_modif,
             filename = paste0(root_for_export,"RiskFunctions/d=",
                length(l_name),"_pts_input_HT_model",
                Z,".png"),
             width=10,height=6)
    }
    
    ### Hull requires a lot of time...
    CHULL<-NA
    RUN_chull<-Params_risk_Function[["RUN_Convex_Hull"]]
    if(RUN_chull){
      GG_hull<-Analysis_diag_HTawn_Conv_Hull(
        Theta_opt = Theta_opt,
        vect_name_variables = l_variables)+
        Common_theme+theme(axis.text=element_text(size=12))
      ggsave(plot=GG_hull,
             filename = paste0(root_for_export,
                               "RiskFunctions/d=",length(l_variables),
                               "_Convex_huls_Htawn.png"),
             width=8,height=4)
    }
    
    prop_cl<-Q_from_mixt
    if(d==3){
      W_scenarios<-rep(Q_thresh,d)
      sub_vector<-apply(Vect_l_function[,c(1:2)],
                        MARGIN = 1,FUN = max)
      new_weight<-mean(as.numeric(sub_vector>Threshold
                              ))-Q_from_mixt[IREF]
      W_scenarios[IREF]<-Q_from_mixt[IREF]
      W_scenarios[2]<-new_weight
      W_scenarios[3]<-(1-Q_thresh)-(new_weight+Q_from_mixt[IREF])
    }
    prop_cl<-W_scenarios
  }
  
  MATRICE_SCORES<-c()
  l_fonction<-list()
  l_Mat_moyenne<-list()
  l_sigma<-list()
  list_SHAPE_OBS<-list()
  d<-ncol(Vect_l_function)
  
  L<-length(l_variables)-1
  for(NameVAR in l_variables[1:L]){
    Name_for_export<-paste0(root_for_export,NameVAR,"_")
  }
  Name_for_export<-paste0(Name_for_export,"_",l_variables[length(l_variables)])
  ## (2) Approximate the distribution of angular tail components --------
  ####################
  
  if(One_PCA_base){
    print("One PCA basis")
    LIST_Mod<-Approach_Angle_One_PCA(LIST_all = list_all,
                                     Name_for_export =root_for_export,
                                     NbScores_Omega = NbScores_Omega,
                                     l_variables = l_variables,
                                     Indices_exts = Indices_exts,
                                     d = d,f_transf = f_transf)
    Inertia_EXPLAINED<-LIST_Mod[["inertia_explained"]]
    print(paste0("inertia explained ", Inertia_EXPLAINED))
  }else{
    print("Several PCA basis")
    LIST_Mod<-Approach_Angle_Mult_PCA(Indices_exts = Indices_exts,
                                      LIST_all = list_all,
                                      root_export = root_for_export,
                                      l_variables = l_variables,
                                      list_nb_scores=list_nb_scores,
                                      f_transf = f_transf)
    NbScores_Omega<-LIST_Mod[["Nb_scores_omega"]]
  }
  list_Frechet_OBS<-LIST_Mod[["L_Frechet"]]
  Scores<-LIST_Mod[["Scores"]]
  Length_T<-LIST_Mod[["Length_TS"]]
  EIGEN_functions<-LIST_Mod[["Eigen_functions"]]
  Unif_coord<-VineCopula::pobs(Scores)
  M_Theta<-Params_risk_Function[["Method_Angle"]]
  if(M_Theta=="VineCop"){
    Matrix_C<-function_Structure_Matrice(NB_dim = NbScores_Omega)
    Model_coords<-VineCopula::RVineCopSelect(Unif_coord,
                                              Matrix = Matrix_C)
  }
  if(M_Theta=="BetaCop"){
    Model_coords<-copula::empCopula(X = Unif_coord,
                      smoothing = "beta")
  }
  if(M_Theta=="GaussMixture"){
    Model_coords<-mclust::densityMclust(Scores,
                                    plot=FALSE)
  }
  Bool_filled<-FALSE
  N_sim<-M
  Coords_simul<-matrix(NA,nrow = M,ncol = NbScores_Omega)
  ### (3) Simulation step, sample from the probabilistic model---------------------
  #################
  
  while(Bool_filled==FALSE){
    
    if(M_Theta!="GaussMixture"){
      if(M_Theta=="VineCop"){
        # Simulations of coordinates  ------------------------------------------------------------
        Simulations_coord<-VineCopula::RVineSim(N=N_sim,
                                                RVM = Model_coords)
      }
      if(M_Theta=="BetaCop"){
        Simulations_coord<-copula::rCopula(N_sim, 
                                           copula=Model_coords)
      }
      
      colnames(Simulations_coord)<-1:ncol(Simulations_coord)
      Coords_ech_orig<-matrix(NA,nrow=N_sim,
                              ncol=NbScores_Omega)
      for(j in c(1:ncol(Simulations_coord))){
        sortie<-quantile(Scores[,j],Simulations_coord[,j])
        Coords_ech_orig[,j]<-as.numeric(sortie)
      }
    }
    if(M_Theta=="GaussMixture"){
      Coords_ech_orig<-Simul_from_MixtureGauss(n = N_sim,
                    Object_dens_clust = Model_coords,
                    NB_dim_PCA = ncol(Scores))
    }
    if(One_PCA_base){
      Shape_Omega_simul<-Simul_Omega_One_PCA_base(list_Mod_One_PCA = LIST_Mod,
                                                  Simul_coords =Coords_ech_orig 
                                                  ,NbScores_Omega =NbScores_Omega,
                                                  f_transf_inv = f_transf_inv)
    }else{
      Shape_Omega_simul<-Simul_Omega_Mult_PCA_base(list_Mod_Mult_PCA =LIST_Mod,
                                                   M = N_sim,d = d,list_nb_scores = list_nb_scores,
                                                   Simul_coords =Coords_ech_orig,
                                                   f_transf_inv = f_transf_inv)
    }
    
    ##### Compar simul of theta with the reality
    scale_frechet_d<-Params_risk_Function[["scale_frechet_d"]]
    if(Approach=="AD"){
      if(Params_risk_Function[["parametric"]]==FALSE){
        UNIF<-exp(-Vect_l_transf^(-1))
        EXP<--log(1-UNIF)
        Sim_l<-Simul_RF_max_NP(Params_risk_Function = Params_risk_Function,
                               M=N_sim,QTH = Q_thresh,
                               Vect_l_function = EXP)
        Lg<-apply(X =EXP,MARGIN = 1,FUN = max)
        Threshold_exp<-quantile(Lg,Q_thresh)
        Sim_l_tf<-pexp(Sim_l+Threshold_exp)
      }else{ 
        Sim_l<-Simulation_gParetoP(Params_risk_Function=Params_risk_Function,
                                   d=ncol(Vect_l_function),
                                   M=N_sim,theta_opt=Theta_opt)
        ### From Legrand
        ### Warning) definition by continuity if null values are given
        Sim_l<-Sim_l*Seuil_Margins_Pareto
        Sim_l_tf<-apply(Sim_l,MARGIN=2,FUN = evd::pgpd,loc=0,
              scale=1, shape=1)
      }
    }
   
    if(Approach=="HTawn"){
      # Rejection sampling approach
      Sim_l_tf<-Simul_from_Htawn(model_mex_all = Model_Htawn_all,
               n_sim = N_sim,d = d,Name_vars = l_variables,
               prop_class = prop_cl,ind_ref=IREF,
               thresh_censor = Threshold,
               Q_sim = DQU_each_scenario)
    }
    Sim_l<-Sim_l_tf
    ### Extreme individuals in the original scale Pareto(shape param)
    Excedents_lprime_orig<-Vect_l_function[Indices_exts,]
    ### Chimeas for sim L (vs) obs L
    png(filename= paste0(root_for_export,"RisKfunctions/chi_meas_sim_obs",
                         Params_risk_Function[["general_option"]],".png"),
        width=750,250)
    par(mfrow=c(1,2))
    POT::chimeas(Excedents_lprime_orig[,c(1:2)],which=1)
    POT::chimeas(Sim_l[,c(1:2)],which=2)
    par(mfrow=c(1,1))
    dev.off()
    ####
    Sim_l<-as.data.frame(Sim_l)
    colnames(Sim_l)<-l_variables
    if(length(l_variables)>2){
      list_R_GM<-list()
      Avail_combinations<-t(utils::combn(x = length(l_variables)
                                         ,m = 2))
      for(j in c(1:nrow(Avail_combinations))){
        PAIR<-Avail_combinations[j,]
        NAMES_avail<-l_variables[PAIR]
        Exces_sub<-Excedents_lprime_orig[,PAIR]
        sim_lsub<-Sim_l[,PAIR]
        Df_combs_sim_exts<-rbind.data.frame(Exces_sub,
                                            sim_lsub)
        
        colnames(Df_combs_sim_exts)<-c("V1","V2")
        Df_combs_sim_exts$Legend<-c(rep("data",
                                        nrow(Exces_sub)),
                                    rep("simulations",
                                        nrow(sim_lsub)))
        Z_gg<-PAIR[1]
        Z_gg2<-PAIR[2]
        XLAB<-expression(R[M~","~Z_gg]^g)
        XLAB<-as.expression(do.call('substitute', list( XLAB[[1]], 
                                                        list(Z_gg=Z_gg))))
        YLAB<-expression(R[M~","~Z_gg2]^g)
        YLAB<-as.expression(do.call('substitute', 
                                    list( YLAB[[1]], list(Z_gg2=Z_gg2))))
        GG_RG_sim_vs_obs<-ggplot(data=Df_combs_sim_exts,
                                 aes(x=V1,y=V2))+
          geom_point(aes(shape=Legend,col=Legend,
                         size=Legend))+
          geom_xsidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                            alpha = 0.5)+
          geom_ysidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                            alpha = 0.5)+
          scale_y_continuous(transform = "log10")+
          scale_x_continuous(transform = "log10")+
          xlab(XLAB)+
          ylab(YLAB)+
          scale_shape_manual(values = c("simulations"=17,
                                        "data"=19))+
          scale_size_manual(values=c("simulations"=0.75,
                                     "data"=1.5))+
          scale_color_manual(values=cols_ggplot)+
          scale_fill_manual(values=cols_ggplot)+
          guides(fill="none")
        if(j==1){
          GG_root<-GG_RG_sim_vs_obs+
            Common_theme+theme(
                legend.direction = "horizontal")
          Combn_RGM<-GG_RG_sim_vs_obs+
            Common_theme+theme(
                  legend.direction = "horizontal",
                  legend.position = "none")
          LGD_RGM<-get_legend(GG_root)
        }else{
          Rest<-j-1
          Transit<-GG_RG_sim_vs_obs+
            Common_theme+
            theme(
                  legend.position = "none",
                  legend.direction = "horizontal")
          Combn_RGM<-plot_grid(Combn_RGM,
                Transit,rel_widths = c(Rest,1))
        }
        
          
      }
    }else{
      Df_combs_sim_exts<-rbind.data.frame(Excedents_lprime_orig,
                                          Sim_l)
      colnames(Df_combs_sim_exts)<-c("V1","V2")
      Df_combs_sim_exts$Legend<-c(rep("data",
                                      nrow(Excedents_lprime_orig)),
                                  rep("simulations",
                                      nrow(Sim_l)))
      Z_gg<-1
      Z_gg2<-2
      XLAB<-expression(R[M~","~Z_gg]^g)
      XLAB<-as.expression(do.call('substitute', list( XLAB[[1]], 
                                                      list(Z_gg=Z_gg))))
      YLAB<-expression(R[M~","~Z_gg2]^g)
      YLAB<-as.expression(do.call('substitute', 
                                  list( YLAB[[1]], list(Z_gg2=Z_gg2))))
      GG_RG_sim_vs_obs<-ggplot(data=Df_combs_sim_exts,
                               aes(x=V1,y=V2))+
        geom_point(aes(shape=Legend,col=Legend,
                       size=Legend))+
        geom_xsidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                          alpha = 0.5)+
        geom_ysidedensity(data=Df_combs_sim_exts,aes(fill=Legend), 
                          alpha = 0.5)+
        scale_y_continuous(transform = "log10")+
        scale_x_continuous(transform = "log10")+
        xlab(XLAB)+
        ylab(YLAB)+
        scale_shape_manual(values = c("simulations"=17,
                                      "data"=19))+
        scale_size_manual(values=c("simulations"=0.75,
                                   "data"=1.5))+
        scale_color_manual(values=cols_ggplot)+
        scale_fill_manual(values=cols_ggplot)+
        guides(fill="none")+
        Common_theme
      Whole_SIM_L<-Df_combs_sim_exts
      ggsave(filename= paste0(root_for_export,"RisKfunctions/Rg_sim_vs_obs_",
                              Params_risk_Function[["general_option"]],
                              "_d=",length(l_name),".png"),
             plot=Whole_SIM_L,
             width=8,height=6)
    }
    ### Export simulated radial vectors -------------
    ###############
    LG<-plot_grid(LGD_RGM)
    Combn_leg<-plot_grid(Combn_RGM,
                         LG,nrow=2,
                         ncol=1,rel_heights =c(10,1))
    
    Combn_leg<-ggdraw()+
      draw_plot(Combn_leg)+
      theme(plot.background = element_rect(fill = "white", 
                                           color = NA))
    File_rmg<-paste0(root_for_export,"RisKfunctions/Rg_sim_vs_obs_",
                     Params_risk_Function[["general_option"]],
                     "_d=",length(l_name),".png")
    ggsave(plot = Combn_leg,filename = File_rmg,
      height=6,width=10)

    list_shapes<-list()
    list_candidats<-list()
    DF_simul_lprime<-c()
    NB_Forcg_cond<-length(l_variables)
    Beg<-0
    list_ind_pos<-list()
    for(k in c(1:NB_Forcg_cond)){
      name_s<-l_variables[k]
      L_prime_simul<-Sim_l[,k]
      #Deconcatenate the Omega----------------
      #####
      Beg<-Beg+1
      END<-Beg+Length_T[[name_s]]-1
      Shape_forcing<-Shape_Omega_simul[,Beg:END]
      L2_shape_name<-apply(Shape_forcing,MARGIN = 1,
                           FUN = calcul_norm_L2)
      Shape_forcing_std<-t(t(Shape_forcing)%*%diag(L2_shape_name^(-1)))
      list_shapes[[name_s]]<-Shape_forcing_std
      Z_varj<-t(t(Shape_forcing_std)%*%(diag(L_prime_simul)))
      list_candidats[[name_s]]<-Z_varj
      DF_simul_lprime<-c(DF_simul_lprime,
            apply(X = Z_varj,MARGIN = 1,
                    FUN = calcul_norm_L2))

      indicatrice_pos<-which(apply(X=Z_varj,
                               FUN=fonction_trajectoire_positive,
                               MARGIN = 1)==TRUE)
      Zsub<-Z_varj[indicatrice_pos,]
      Unif<-apply(Z_varj[indicatrice_pos,],
                  MARGIN=2,FUN = evd::pgpd,loc=0,
                      scale=Shape_evd, shape=Shape_evd)
      PTO<-(1-Unif)^(-1)
      seuil_marg<-1+10^(-5)
      
      #### select the simulated multivariate time series above
      #### the threshold of the GPD threshold
      
      indicatrice_pos2<-which(apply(X=(PTO-seuil_marg),
                FUN=fonction_trajectoire_positive,
                MARGIN = 1)==TRUE)
      indicatrice_pos<-indicatrice_pos[indicatrice_pos2]
      if(k==1){
        ### initialisation
        list_ind_pos[[k]]<-indicatrice_pos
      }else{
        ### intersection of intersection etc...
        list_ind_pos[[k]]<-intersect(list_ind_pos[[k-1]],
                                     indicatrice_pos)
      }
      
      Beg<-END
    }
    ref<-Array_simul[1,,1]
    Nb_filled<-length(which(!is.na(ref)==TRUE))
    Inds_final_chosen<-list_ind_pos[[NB_Forcg_cond]]
    Nb_final_chosen<-length(Inds_final_chosen)
    Beg_filled<-Nb_filled+1
    Gap<-M-(Nb_filled+Nb_final_chosen)
    End_filled<-min(Nb_filled+Nb_final_chosen,M)
    Lag<-End_filled-Beg_filled+1
    for(k in c(1:length(l_variables))){
      name_variable<-l_variables[k]
      Z_vark<-list_candidats[[name_variable]][Inds_final_chosen,]
      Sub_coords<-Coords_ech_orig[Inds_final_chosen,]
      Sub_coords<-Sub_coords[c(1:Lag),]

      UNIF_k<-apply(Z_vark,MARGIN=2,
                  FUN = evd::pgpd,loc=0,
                  scale=Shape_evd, shape=Shape_evd)
      UNIF_k<-UNIF_k[c(1:Lag),]
      ### Unif scale
      D_full<-dim(Array_simul)[3]
      D<-ncol(UNIF_k)
      Filled_indexes<-c(Beg_filled:End_filled)
      Nfound<-length(Filled_indexes)
      ### Fill with 0 if nb_col<full case
      Array_simul[k,Filled_indexes,c(1:D_full)]<-matrix(0,
                              nrow = length(Nfound),
                              ncol=D_full)
      Array_simul[k,Filled_indexes,c(1:D)]<-UNIF_k
      Coords_simul[Filled_indexes,]<-Sub_coords
    }
    ### See if the object is completed. 
    new_ref<-Array_simul[1,,1]
    Nb_filled_update<-length(which(!is.na(new_ref)==TRUE))
    ### if True, stop the loop
    if(Nb_filled_update==M){
      Bool_filled<-TRUE
    }
  }

  # Conversion of the simulated multivariate TS --------------------------------------------------------------
  # ------------------------------------------------------------------------
  list_Frechet_SIMUL<-list()
  list_simul<-list()
  list_obs_exts<-list()
  print(Length_T)
  
  # Since two methods are possible to reconvert the MV TS, 
  # we look at the one used in the code
  
  if(result_transformation[["type_transfo"]]=="_Mixt_transf_"){
    for(k in c(1:length(l_variables))){
      name_variable_for_conv<-l_variables[k]
      D<-Length_T[[name_variable_for_conv]]
      Simul_whole_kunif<-Array_simul[k,,c(1:D)]
      
      ### Unif-->Pareto to compare in the scale of obs.
      # If GPD(1,0,1) in conversion
      Conv_for_analyse<-apply(Simul_whole_kunif,MARGIN=2,
            FUN = evd::qgpd,loc=0,
            scale=Shape_evd, 
            shape=Shape_evd)
      Times_used<-Times_per_FCG[[name_variable_for_conv]]
      DF<-as.data.frame(Conv_for_analyse)
      colnames(DF)<-Times_used
      list_Frechet_SIMUL[[name_variable_for_conv]]<-DF
      INDICES_ACP<-1:nrow(Z_varj)
      # Reconversion in the correct scale--------------------------------------
      Theta_EXTGPD_k<-list_all[[name_variable_for_conv]]$params_transfo
      K<-Theta_EXTGPD_k$K
      LEVT<-Theta_EXTGPD_k$LEVT

      DFF<-as.data.frame(Simul_whole_kunif)
      colnames(DFF)<-Times_used
      Variables_reconversion_ACP<-lapply(Times_used,
                    function_reconversion_Pareto_fromDframes,
                     Df_K=K,variable_uplift=DFF,
                     Df_evt=LEVT)
      Variables_reconversion_ACP<-do.call(cbind.data.frame,
                     Variables_reconversion_ACP)
      colnames(Variables_reconversion_ACP)<-Times_used
      list_simul[[name_variable_for_conv]]<-Variables_reconversion_ACP
      
      # Extreme observations -----------------------------------------------------
      obs_ext<-result_transformation$orig[[name_variable_for_conv]][Indices_exts,]
      colnames(obs_ext)<-Times_used
      list_obs_exts[[name_variable_for_conv]]<-obs_ext
    }
  }else{
    for(k in c(1:length(l_variables))){
      name_variable_for_conv<-l_variables[k]
      D<-Length_T[[name_variable_for_conv]]
      
      Simul_whole_kunif<-Array_simul[k,,c(1:D)]
      #Conv_for_analyse<-(-1)*(log(Simul_whole_kunif))^(-1)
      ### Unif-->Pareto to compare in the scale of obs.
      Conv_for_analyse<-apply(Simul_whole_kunif,MARGIN=2,
                              FUN = evd::qgpd,loc=0,
                              scale=Shape_evd, 
                              shape=Shape_evd)
      Times_used<-Times_per_FCG[[name_variable_for_conv]]
      INDICES_ACP<-1:nrow(Simul_whole_kunif)
      ### Unif-->Frechet to compare in Frechet scale of obs.
      DF<-as.data.frame(Conv_for_analyse)
      colnames(DF)<-Times_used
      list_Frechet_SIMUL[[name_variable_for_conv]]<-DF
      NO_ACP<-lapply(INDICES_ACP,FUN = fnct_select_colonne,
                     df=Simul_whole_kunif)
      
      # Reconversion in the correct scale--------------------------------------
      ### Change the colnames using correct times
      Theta_EXTGPD_k<-list_all[[name_variable_for_conv]]$params_transfo
      Variables_reconversion_ACP<-lapply(NO_ACP,
                    function_reconversion_Unif_EXTGPD,
                    Theta_k=Theta_EXTGPD_k)
      
      Variables_reconversion_ACP<-do.call(rbind,
                      Variables_reconversion_ACP)
      colnames(Variables_reconversion_ACP)<-Times_used
      list_simul[[name_variable_for_conv]]<-Variables_reconversion_ACP
      
      # Extreme observations -----------------------------------------------------
      obs_ext<-result_transformation$orig[[name_variable_for_conv]][Indices_exts,]
      colnames(obs_ext)<-Times_used
      list_obs_exts[[name_variable_for_conv]]<-obs_ext
    }
    
  }
  if(M_Theta=="GaussMixture"){
    Object_modelisation<-Model_coords
  }
  else{
    Object_modelisation<-Model_coords
  }
  return(list("Param_found"=Theta_opt,"simul"=list_simul,
              "obs_exts"=list_obs_exts,"Indices_exts"=Indices_exts,
              "coords_simul"=Coords_simul, "coords_data"=Scores,
              "Frechet_normal"=list("obs"=list_Frechet_OBS,
                                    "simul"=list_Frechet_SIMUL),
              "EIGEN_functions"= EIGEN_functions,
              "threshold"=Threshold,
              "Vect_RiskF_transf"=Vect_l_transf,
              "model_EGPD"=L_Transf_vect_l,
              "mqu"=MQU,"vect_transf_EGPD"=Vect_Pareto_margins,
              "Model_Coords"=Model_coords,
              "CHULL_HTAWN"=CHULL))
}

#' function_reconversion_Unif_EXTGPD
#'
#' @param Simul_Unif vector. Simulated vector with uniform margins.
#' @param Theta_k list. List of EGPD estimators found at each time t.
#'
#' @return Simulated vector with EGPD margins.
#' @export
#'
#' @examples
function_reconversion_Unif_EXTGPD<-function(Simul_Unif,Theta_k){
  
  L_dim<-length(Simul_Unif)
  individu_ext<-sapply(c(1:L_dim),function_reconv_each_dim_EGPD,Simul_Unif=Simul_Unif,
                       Theta_k=Theta_k)
  return(individu_ext)
}
#' function_reconv_each_dim_EGPD
#'
#' @param Simul_Unif vector. Simulated vector with uniform margins.
#' @param Theta_k list. List of EGPD estimators found at each time t.
#' @param index_dim int. Time index. 
#'
#' @return Value at time t of the vector at EGPD scale. 
#' @export
#'
#' @examples
function_reconv_each_dim_EGPD<-function(Simul_Unif,Theta_k,index_dim){
  
  value_unif_x<-Simul_Unif[index_dim]
  Theta_kt<-Theta_k[[index_dim]]
  Nu_t<-Theta_kt[["nu"]]
  Gamma_t<-Theta_kt[["mu"]]
  Sigma_t<-Theta_kt[["sigma"]]
  #Apply inverse of G
  value_G_x<-value_unif_x^(1/Nu_t)
  #Apply inverse of GPD cd
  value_EGPD_t<-((1-value_G_x)^(-Gamma_t)-1)*(Sigma_t/Gamma_t)
  return(value_EGPD_t)
}


Function_one_couple_l<-function(Params_Biv,z,Th,one_u){

  realisation<-Function_mixture(Params_Biv =Params_Biv ,u = one_u,
                   Th =Th,z = z)
  realisation_y<-Simul_cond_one(input_x = realisation,
                 Params_biv = Params_Biv,
                 Th = Th,
                 z = z)
  return(c(realisation,realisation_y))
}



#' Couple_s_t_chi_measure
#'
#' @param t int. First time of the pair.  
#' @param s int. Second time of the pair. 
#' @param matrix_ matrix. Observations matrix.
#' @param tail_quantile float. High threshold in the uniform scale. 
#'
#' @return Extremal correlation coefficient Chi for the pair at (t,s) 
#' @export
#'
#' @examples
Couple_s_t_chi_measure<-function(t,s,matrix_,tail_quantile){
  
    Sub_matrix<-matrix_[,c(t,s)]
    Min_st<-apply(X = Sub_matrix,MARGIN = 1,FUN = min)
    prob_intersection<-mean(as.numeric(Min_st>tail_quantile))
    return(prob_intersection/(1-tail_quantile))
}
#'Matrix_chi_measure
#'
#' @param Data_unif matrix. Observation with uniform margins. 
#' @param tail_quantile float. High threshold in the uniform scale. 
#'
#' @return Matrix of extremal correlation coefficient Chi per pair (t,s) 
#' @export
#'
#' @examples
Matrix_chi_measure<-function(Data_unif,tail_quantile){
  
  matrix_chi<-matrix(0,nrow = ncol(Data_unif),
                     ncol=ncol(Data_unif))
  for(i in c(1:ncol(matrix_chi))){
    sup_ij<-c(i:d)
    Values_chi<-sapply(X = sup_ij,FUN = Couple_s_t_chi_measure,
           t=i,matrix_=Data_unif,tail_quantile=tail_quantile)
    matrix_chi[i,sup_ij]<-Values_chi
  }
  Comb<-t(matrix_chi)+matrix_chi
  diag(Comb)<-diag(Comb)/2
  return(Comb)
}
ALD_reg_XY<-function(data_XY,basis_functionGAM,
                     Level_Thresh,Sim_X){
  colnames(data_XY)<-c("X_regressor","Y_target")
  mu_data<-colMeans(data_XY)
  sd_data<-apply(X = data_XY,MARGIN = 2,
                 FUN = sd)
  if(is.null(dim(Sim_X))){
    N_sim<-length(Sim_X)
  }else{
    N_sim<-nrow(Sim_X)
  }
  
  if(is.null(dim(Sim_X))){
    Sim_X1<-scale(Sim_X,center = mu_data[1],
                  scale = sd_data[1])
    Sim_X2<-scale(Sim_X,center = mu_data[2],
                  scale = sd_data[2])
    N_sim<-length(Sim_X)
    Sim_for_output<-matrix(
      rep(Sim_X,2),nrow = N_sim,byrow = FALSE)
  }else{
    Sim_X1<-scale(Sim_X[,1],center = mu_data[1],
                  scale = sd_data[1])
    Sim_X2<-scale(Sim_X[,2],center = mu_data[2],
                  scale = sd_data[2])
    N_sim<-nrow(Sim_X)
    Sim_for_output<-Sim_X
  }
  Sim_X1<-Sim_for_output[,1]
  Sim_X2<-Sim_for_output[,2]
  tau_used<-rep(Level_Thresh,2*N_sim)
  if(basis_functionGAM=="cr"){
    fmla_ald<-paste0('Y_target ~ s(X_regressor,bs="cr" )')
  }
  if(basis_functionGAM=="tp"){
    fmla_ald<-paste0('Y_target ~ s(X_regressor,bs="tp" )')
  }
  
  ###Normalising first coordinate
  data_1<-data_XY
  data_1[,1]<-scale(data_1[,1])
  Ald_model<-evgam::evgam(as.formula(fmla_ald), 
              data=data_1, 
              family="ald", 
              ald.args=list(tau=Level_Thresh))

  ### Inverse 
  data_flip<-data_XY[,c(2,1)]
  data_flip[,1]<-scale(data_flip[,1])
  colnames(data_flip)<-colnames(data_1)
  Ald_model_flip<-evgam::evgam(as.formula(fmla_ald), 
              data=data_flip, 
              family="ald", 
              ald.args=list(tau=Level_Thresh))
  ### Predictions from simulations (first)
  ##########
  Y_pred_Tau_withse<-predict(Ald_model, 
               newdata=list(X_regressor=Sim_X1),
               type="response",se.fit=TRUE)
  ### Use confidence band
  Y_pred_Tau<-Y_pred_Tau_withse[["fitted"]][["location"]]
  Sd_found<-Y_pred_Tau_withse[["se.fit"]][["location"]]
  Bounds<-t(sapply(c(1:length(Y_pred_Tau)),
         function(x,vect_reg){
      MU<-vect_reg[x]
      Sd<-Sd_found[x]
      Q_min_max<-qnorm(p = c(0.025,0.975),mean = MU,
                sd = Sd)
      return(Q_min_max)
         },vect_reg=Y_pred_Tau))
  data_1case<-cbind(Sim_for_output[,1],
                    Y_pred_Tau,Bounds)
  ### Predictions from simulations (second)
  ##########
  Y_pred_Tau_flip_withse<-predict(Ald_model_flip, 
                newdata=list(X_regressor=Sim_X2),
                type="response",
                se.fit=TRUE)
  ### Use confidence band
  Y_pred_Tau_flip<-Y_pred_Tau_flip_withse[["fitted"]][["location"]]
  Sd_found_flip<-Y_pred_Tau_flip_withse[["se.fit"]][["location"]]
  Bounds_flip<-t(sapply(c(1:length(Y_pred_Tau_flip)),
       function(x,vect_reg){
         MU<-vect_reg[x]
         Sd<-Sd_found_flip[x]
         Q_min_max<-qnorm(p = c(0.025,0.975),mean = MU,
                          sd = Sd)
         return(Q_min_max)},
        vect_reg=Y_pred_Tau_flip)
  )
  
  data_2case<-cbind(Y_pred_Tau_flip,
                    Sim_for_output[,2],Bounds_flip)
  Descrip_case<-rep("first",N_sim)
  Descrip_case_flip<-rep("second",N_sim)
  Descrip_case_wh<-c(Descrip_case,
                Descrip_case_flip)
  Descrip_case_wh<-cbind(tau_used,Descrip_case_wh)
  Couple_XY<-rbind(data_1case,
                   data_2case)
  Df_obtained_full<-cbind.data.frame(Couple_XY,
                    Descrip_case_wh)
  colnames(Df_obtained_full)<-c("X_reg","Y_target",
                                "ymin","ymax","Tau",
                                "origin")
  return(Df_obtained_full)

}
craft_mex_tform<-function (x, margins, r = NULL, method = "mixture", divisor = "n+1", 
                           na.rm = TRUE) 
{
  margins <- list(casefold(margins), p2q = switch(casefold(margins), 
                                                  gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
                                                                                                                     0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
                                                                                                                                                                        gumbel = function(q) exp(-exp(-q)), 
                                                                                                                                                                        laplace = function(q) ifelse(q < 
                                                                                                                                                                                                       0, exp(q)/2, 1 - 0.5 * exp(-q))))
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
laplace_function<-function(p) {ifelse(p < 
                                        0.5, log(2 * p), 
                              -log(2 * (1 - p)))
}

### Plot functions density Mclust

#' CRAFt_plot_density_D
#'
#' @param x 
#' @param data 
#' @param nlevels 
#' @param levels 
#' @param prob 
#' @param points.pch 
#' @param points.col 
#' @param points.cex 
#' @param gap 
#' @param cex.annot 
#' @param ... 
#'
#' @return Result of the slightly modified 
#' mclust::plot_density_D function. 
#' @export
#'
#' @examples
CRAFt_plot_density_D<-function (x, data = NULL, nlevels = 11, levels = NULL, prob = c(0.25, 
                                                                                      0.5, 0.75), points.pch = 1, points.col = 1, points.cex = 0.8, 
                                gap = 0.2, cex.annot,...) 
{
  
  object <- x
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$points.pch <- mc$points.col <- mc$points.cex <- mc$gap <- NULL
  mc$nlevels <- nlevels
  mc$levels <- levels
  mc$prob <- prob
  if (!is.null(mc$type)) 
    if (mc$type == "level") 
      mc$type <- "hdr"
  if (is.null(data)) {
    data <- mc$data <- object$data
    addPoints <- FALSE
  }
  else {
    data <- as.matrix(data)
    stopifnot(ncol(data) == ncol(object$data))
    addPoints <- TRUE
  }
  nc <- object$d
  oldpar <- par(mfrow = c(nc, nc), mar = rep(gap/2, 4), oma = rep(3, 
                                                                  4), no.readonly = TRUE)
  on.exit(par(oldpar))
  for (i in seq(nc)) {
    for (j in seq(nc)) {
      if (i == j) {
        plot(data[, c(i, j)], type = "n", xlab = "", 
             ylab = "", axes = FALSE)
        text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
             colnames(data)[i], cex =cex.annot, adj = 0.5)
        box()
      }
      else {
        par <- object$parameters
        if (is.null(par$pro)) 
          par$pro <- 1
        par$mean <- par$mean[c(j, i), , drop = FALSE]
        par$variance$d <- 2
        sigma <- array(dim = c(2, 2, par$variance$G))
        for (g in seq(par$variance$G)) sigma[, , g] <- par$variance$sigma[c(j, 
                                                                            i), c(j, i), g]
        par$variance$sigma <- sigma
        par$variance$Sigma <- NULL
        par$variance$cholSigma <- NULL
        par$variance$cholsigma <- NULL
        mc$parameters <- par
        mc$data <- object$data[, c(j, i)]
        mc$axes <- FALSE
        mc[[1]] <- as.name("surfacePlot")
        eval(mc, parent.frame())
        box()
        if (addPoints & (j > i)) 
          points(data[, c(j, i)], pch = points.pch, col = points.col, 
                 cex = points.cex)
      }
      if (i == 1 && (!(j%%2))) 
        axis(3)
      if (i == nc && (j%%2)) 
        axis(1)
      if (j == 1 && (!(i%%2))) 
        axis(2)
      if (j == nc && (i%%2)) 
        axis(4)
    }
  }
  invisible()
}
#' Title
#'
#' @param x 
#' @param G 
#' @param modelNames 
#' @param symbols 
#' @param colors 
#' @param xlab 
#' @param ylab 
#' @param legendArgs 
#' @param cex.lab.plot 
#' @param cex.axis.plot 
#' @param ... 
#'
#' @return Result of the slightly modified 
#' mclust::plot_BIC function. 
#' @export
#'
#' @examples
CRAFt_plot_BIC<-function (x, G = NULL, modelNames = NULL, symbols = NULL, colors = NULL, 
                          xlab = NULL, ylab = "BIC", legendArgs = list(x = "bottomright", 
                                                                       ncol = 2, cex = 1, inset = 0.01), 
                          cex.lab.plot,cex.axis.plot,...) 
{
  
  args <- list(...)
  if (is.null(xlab)) 
    xlab <- "Number of components"
  subset <- !is.null(attr(x, "initialization")$subset)
  noise <- !is.null(attr(x, "initialization")$noise)
  ret <- attr(x, "returnCodes") == -3
  legendArgsDefault <- eval(formals(plot.mclustBIC)$legendArgs)
  legendArgs <- append(as.list(legendArgs), legendArgsDefault)
  legendArgs <- legendArgs[!duplicated(names(legendArgs))]
  n <- ncol(x)
  dnx <- dimnames(x)
  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  if (is.null(modelNames)) 
    modelNames <- dimnames(x)[[2]]
  if (is.null(G)) 
    G <- as.numeric(dimnames(x)[[1]])
  if (is.null(symbols)) {
    colNames <- dimnames(x)[[2]]
    m <- length(modelNames)
    if (is.null(colNames)) {
      symbols <- if (m > 9) 
        LETTERS[1:m]
      else as.character(1:m)
      names(symbols) <- modelNames
    }
    else {
      symbols <- mclust.options("bicPlotSymbols")[modelNames]
    }
  }
  if (is.null(colors)) {
    colNames <- dimnames(x)[[2]]
    if (is.null(colNames)) {
      colors <- 1:m
      names(colors) <- modelNames
    }
    else {
      colors <- mclust.options("bicPlotColors")
      if (!is.null(names(colors)) & !any(names(colors) == 
                                         "")) 
        colors <- colors[modelNames]
    }
  }
  x <- x[, modelNames, drop = FALSE]
  ylim <- if (is.null(args$ylim)) 
    range(as.vector(x[!is.na(x)]))
  else args$ylim
  matplot(as.numeric(dnx[[1]]), x, type = "b", xaxt = "n", 
          xlim = range(G), ylim = ylim, pch = symbols, col = colors, 
          lty = 1, xlab = xlab, ylab = ylab, main = "",
          cex.lab=cex.lab.plot,cex.axis=cex.axis.plot)
  axis(side = 1, at = as.numeric(dnx[[1]]),
       cex.axis=cex.axis.plot)
  if (!is.null(legendArgs)) {
    do.call("legend", c(list(legend = modelNames, col = colors, 
                             pch = symbols), legendArgs))
  }
  invisible(symbols)
}
Launch_test<-function(BIC,cex.lab.plot,legendArgs,cex.axis.plot,...){
  Sub_f<-function(CEX,LGD,CEX.AXIS,...){
    CRAFt_plot_BIC(x = BIC,cex.lab.plot = CEX,
                   legendArgs = LGD,cex.axis.plot = CEX.AXIS,...)
  }
  return(Sub_f(CEX=cex.lab.plot,LGD=legendArgs,
               CEX.AXIS=cex.axis.plot))
}
