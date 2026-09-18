
#' compute the extremal correlation coefficient between 
#' two observation times (t-th and s-th) of an univariate time series. 
#'
#' @param Matrix_df_unif Matrix. Univariate time series with uniform margins.  
#' @param Order_quantile Float. Threshold used to compute the extremal correlation coefficient. 
#' @param t first time step.
#' @param s second time steps
#'
#' @return
#' @export
#'
#' @examples
Chi_measure_analysis_pair<-function(Matrix_df_unif,Order_quantile,
                                    t,s){
  
  Couple<-Matrix_df_unif[,c(t,s)]
  Max_<-apply(X = Couple,MARGIN = 1,FUN = max)
  prob_max<-mean(as.numeric(Max_<=Order_quantile))
  chi <- 2 - log(prob_max)/log(Order_quantile)
  return(chi)
}
#' Compute the asymptotic dependence test 
#' for a given pair of time steps and observation matrix. 
#'
#' @param Matrix_df Matrix. Univariate time series with uniform margins.  
#' @param Mat_cthresh Matrix. Matrix of threshold values to apply the AD test. 
#' @param t first time step.
#' @param s second time steps
#'
#' @return Dataframe with pvalue and exceedance rate. 
#' @export
#'
#' @examples
AD_test_pair<-function(Matrix_df,Mat_cthresh,
                       t,s){
  
  Couple<-Matrix_df[,c(t,s)]
  Cts<-Mat_cthresh[t,s]
  Test_ad<-taildep.test(Couple[,1], Couple[,2],cthresh = Cts)
  PVAL<-as.numeric(Test_ad$p.value)
  Prop<-as.numeric(Test_ad$parameter[4])*(10)^(-2)
  df<-as.data.frame(cbind(PVAL,Prop))
  colnames(df)<-c("pval","prop")
  return(df)
}
#' Compute the asymptotic dependence test 
#' for an pair of time steps and observation matrix. 
#'
#' @param Matrix_df Matrix. Univariate time series with uniform margins.  
#' @param Mat_cthresh Matrix. Matrix of threshold values to apply the AD test. 
#' @param Filename string. Name of the exported graphic. 
#' @param Name_main String (NA by default). Name of the variable 
#' analysed (used for the title of ggplot2 object). 
#' @param New_breaks_labels vector[float]. Plot option
#' to replace the raw x and y values by the variable name. 
#' @param only_pval 
#' @param THEME_ad_ai 
#'
#' @return
#' @export
#'
#' @examples
AD_test_analysis<-function(Matrix_df,Mat_cthresh,
                           Filename=NA,Name_main=NA,
                           New_breaks_labels=NA,
                           only_pval=NA,THEME_ad_ai){
  
  d<-ncol(Matrix_df)
  Mat_pval_found<-matrix(NA,nrow = d,
                         ncol=d)
  Mat_prop_found<-matrix(NA,nrow = d,
                         ncol=d)
  d_minus<-d-1
  for(i in c(1:d_minus)){
    beg<-i+1
    range_i<-c(beg:d)
    Result_line<-lapply(range_i,FUN = AD_test_pair,
                        t=i,Mat_cthresh=Mat_cthresh,
                        Matrix_df=Matrix_df)
    Correspdg_df<-do.call(rbind.data.frame,
                          Result_line)
    Mat_pval_found[i,range_i]<-as.numeric(Correspdg_df[,"pval"])
    Mat_prop_found[i,range_i]<-as.numeric(Correspdg_df[,"prop"])
  }
  Mat_pval_found<-t(Mat_pval_found)
  rownames(Mat_pval_found)<- as.character(c(1:ncol(Mat_pval_found)))
  colnames(Mat_pval_found) <- as.character(c(1:ncol(Mat_pval_found)))
  m_pval<- melt(Mat_pval_found)
  colnames(m_pval) <- c("Row", "Col", "Value")
  m_pval$category<-rep("pvalue",nrow(m_pval))
  
  ### Proportion
  Mat_prop_found<-t(Mat_prop_found)
  rownames(Mat_prop_found)<- as.character(c(1:ncol(Mat_prop_found)))
  colnames(Mat_prop_found) <- as.character(c(1:ncol(Mat_prop_found)))
  m_prop<- melt(Mat_prop_found)
  colnames(m_prop) <- c("Row", "Col", "Value")
  m_prop$category<-rep("exceedance rate",nrow(m_prop))
  ## NA gives error with ==--> use isTrue to prevent this
  Case_1<-isTRUE(only_pval)
  if(Case_1){
    All_AD_AI<-m_pval
  }else{
    ### Concatenation
    All_AD_AI<-rbind.data.frame(m_pval,m_prop)
  }
  ### Title graph
  Name_main_graph<-paste0("AD test of ",Name_main)
  
  d<-ncol(Matrix_df)
  All_AD_AI$round_Value<-round(All_AD_AI$Value,2)
  GG_complete_AD_AI<-ggplot2::ggplot(All_AD_AI, aes(x = Col, y = Row, fill = Value)) 
  if(!Case_1){
    ### if two elements--> facet_wrap
    GG_complete_AD_AI<-GG_complete_AD_AI+
      facet_wrap(~category) 
  }
  GG_complete_AD_AI<-GG_complete_AD_AI+
    geom_tile()
  if(d<5){
    GG_complete_AD_AI<-GG_complete_AD_AI+
      geom_text(aes(label = round_Value), color = "black",
                size = 3)
  }
  GG_complete_AD_AI<-GG_complete_AD_AI+   
    xlab("")+ylab("")+
    scale_fill_viridis(na.value="lightgrey",
                       alpha = 0.8)+
    # labs(title=Name_main_graph)+
    labs(fill="Legend")+
    THEME_ad_ai+theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          legend.direction="horizontal",
          legend.position="bottom")
  if(sum(!is.na(New_breaks_labels))==length(New_breaks_labels)){
    pb <- ggplot_build(GG_complete_AD_AI)
    y_breaks <-pb$layout$panel_params[[1]]$y$breaks
    New_breaks_labels<-New_breaks_labels[y_breaks]
    GG_complete_AD_AI<-GG_complete_AD_AI+
      scale_x_continuous(breaks = y_breaks,
                         labels = New_breaks_labels)+
      
      scale_y_continuous(breaks = y_breaks,
                         labels = New_breaks_labels)
  }
  ggsave(filename = Filename,
         plot = GG_complete_AD_AI,
         width=12,height = 6)
  return(GG_complete_AD_AI)
}

#' Compute for a given variable the pairwise extremal
#' correlation coefficients. 
#'
#' @param Matrix_df_unif Data matrix with uniform margins. 
#' @param Order_quantile Float. Quantile used for computing 
#' extremal correlation coefficient.
#' @param Filename String (NA by default). Directory path
#' for export. 
#' @param Name_main String (NA by default). Name of the variable 
#' analysed (used for the title of ggplot2 object). 
#'
#' @return Matrix of pairwise correlation coefficient or ggplot2 
#' object or ggplot2 object
#' @export
#'
#' @examples
Chi_measure_analysis<-function(Matrix_df_unif,Order_quantile,
                               Filename=NA,Name_main=NA){
  
  d<-ncol(Matrix_df_unif)
  Mat_Xi_found<-matrix(NA,nrow = d,
                       ncol=d)
  d_minus<-d-1
  for(i in c(1:d_minus)){
    beg<-i+1
    range_i<-c(beg:d)
    Result_line<-sapply(range_i,FUN = Chi_measure_analysis_pair,
                        t=i,Order_quantile=Order_quantile,
                        Matrix_df_unif=Matrix_df_unif)
    Mat_Xi_found[i,range_i]<-Result_line
    
  }
  if(is.na(Filename)){
    return(t(Mat_Xi_found))
  }else{
    mat <-All
    rownames(mat) <- as.character(c(1:ncol(mat)))
    colnames(mat)<- as.character(c(1:ncol(mat)))
    Name_froot<-expression(chi~"matrix of "~ NAME_VAR)
    Name_main_modif<-do.call("substitute", 
                             list(Name_froot[[1]], 
                                  list(NAME_VAR = Name_main)))  
    # Convert matrix to long format
    df <- melt(mat)
    colnames(df) <- c("Row", "Col", "Value")
    GGdefault<-ggplot2::ggplot(df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      xlab("")+ylab("")+
      scale_fill_viridis()+
      labs(title=Name_main_modif)+
      labs(fill="Value")+
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(size=15))
    ggsave(filename = Filename,
           plot = GGdefault,
           width=8,height = 6)
    return(GGdefault)
  }
  
}
#' chi-bar extremal correlation coefficient
#' for a specific pair of components (t,s) of the input matrix
#'
#' @param Matrix_df_unif dataframe[float]. Univariate
#' matrix with uniform margins
#' @param Order_quantile float. Quantile in the uniform scale.
#' @param t int. Coordinate index to analyse. 
#' @param s int. Other coordinate index to analyse. 
#'
#' @return Float.
#' @export
#'
#' @examples
Chi_bar_measure_analysis_pair<-function(Matrix_df_unif,Order_quantile,
                                        t,s){
  
  Couple<-Matrix_df_unif[,c(t,s)]
  Min_<-apply(X = Couple,MARGIN = 1,FUN = min)
  prob_min<-mean(as.numeric(Min_>Order_quantile))
  #From chimeas
  chibar <- 2 * log(1-Order_quantile)/log(prob_min) - 1
  return(chibar)
}

#' Ggplot object displaying the chi-bar extremal correlation coefficient
#' for every pair of components of the input matrix.
#'
#' @param Matrix_df_unif dataframe[float]. Univariate
#' matrix with uniform margins
#' @param Order_quantile float. Quantile in the uniform scale.
#' @param Filename string. Name of the exported graphic. 
#' @param Name_main String (NA by default). Name of the variable 
#' analysed (used for the title of ggplot2 object). 
#'
#' @return Ggplot object. 
#' @export
#'
#' @examples
Chi_bar_measure_analysis<-function(Matrix_df_unif,Order_quantile,
                                   Filename=NA,Name_main=NA){
  
  d<-ncol(Matrix_df_unif)
  Mat_Xibar_found<-matrix(NA,nrow = d,
                          ncol=d)
  d_minus<-d-1
  for(i in c(1:d_minus)){
    beg<-i+1
    range_i<-c(beg:d)
    Result_line<-sapply(range_i,FUN = Chi_bar_measure_analysis_pair,
                        t=i,Order_quantile=Order_quantile,
                        Matrix_df_unif=Matrix_df_unif)
    Mat_Xibar_found[i,range_i]<-Result_line
    
  }
  if(is.na(Filename)){
    return(t(Mat_Xibar_found))
  }else{
    mat <-All
    rownames(mat) <- as.character(c(1:ncol(mat)))
    colnames(mat)<- as.character(c(1:ncol(mat)))
    
    # Convert matrix to long format
    df <- melt(mat)
    input_for_tex<-paste0("$\\bar{\\chi}$ matrix of ",Name_main)
    Name_main_modif<-latex2exp::TeX(input_for_tex)
    colnames(df) <- c("Row", "Col", "Value")
    GGdefault<-ggplot2::ggplot(df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      xlab("")+ylab("")+
      scale_fill_viridis()+
      labs(title=Name_main_modif)+
      labs(fill="Value")+
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(size=15))
    ggsave(filename = Filename,
           plot = GGdefault,
           width=8,height = 6)
    return(GGdefault)
  }
  
}
#' Complete_chi_measure_analysis
#'
#' @param Matrix_df_unif dataframe[float]. Univariate
#' matrix with uniform margins
#' @param Order_quantile float. Quantile in the uniform scale.
#' @param Filename string. Name of the exported graphic. 
#' @param Name_main String (NA by default). Name of the variable 
#' analysed (used for the title of ggplot2 object). 
#' @param New_breaks_labels vector[float]. Plot option
#' to replace the raw x and y values by the variable name. 
#' @param Lim_viridis Limits of the extremal correlation 
#' coefficients in the ggplot objects. 
#'
#' @return
#' @export
#'
#' @examples
Complete_chi_measure_analysis<-function(Matrix_df_unif,Order_quantile,
                                        Filename=NA,Name_main=NA,
                                        New_breaks_labels=NA,
                                        Lim_viridis=NULL){
  
  df_chi<-Chi_measure_analysis(
    Matrix_df_unif = Matrix_df_unif,
    Order_quantile = Order_quantile,Filename = NA,
    Name_main = Name_main)
  rownames(df_chi) <- as.character(c(1:ncol(df_chi)))
  colnames(df_chi) <- as.character(c(1:ncol(df_chi)))
  m_chi<- melt(df_chi)
  colnames(m_chi) <- c("Row", "Col", "Value")
  m_chi$category<-rep("chi",nrow(m_chi))
  
  df_chi_bar<-Chi_bar_measure_analysis(
    Matrix_df_unif = Matrix_df_unif,
    Order_quantile = Order_quantile,Filename = NA,
    Name_main = Name_main)
  colnames(df_chi_bar)<- as.character(c(1:ncol(df_chi_bar)))
  rownames(df_chi_bar)<- as.character(c(1:ncol(df_chi_bar)))
  
  m_chi_bar<- melt(df_chi_bar)
  colnames(m_chi_bar) <- c("Row", "Col", "Value")
  m_chi_bar$category<-rep("chi_bar",
                          nrow(m_chi_bar))
  ### Concatenation
  All_chi<-rbind.data.frame(m_chi,m_chi_bar)
  ### Title graph
  Name_main_graph<-paste0("Correlation matrix of ",Name_main)
  cat_labels<-c("chi"="chi",
                "chi_bar"="bar(chi)")
  d<-ncol(Matrix_df_unif)
  All_chi$round_Value<-round(All_chi$Value,2)
  GG_complete_chi<-ggplot2::ggplot(All_chi, aes(x = Col, y = Row, fill = Value)) +
    facet_wrap(~category,
               labeller=as_labeller(cat_labels,label_parsed))+
    geom_tile()
  if(d<5){
    GG_complete_chi<-GG_complete_chi+
      geom_text(aes(label = round_Value), color = "black",
                size = 3)
  }
  GG_complete_chi<-GG_complete_chi+   
    xlab("")+ylab("")+
    scale_fill_viridis_c(na.value="lightgrey",
                         alpha = 0.8,
                         limits=Lim_viridis)+
    labs(fill="Value")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          strip.text = element_text(size=14),
          axis.text  = element_text(size=12),
          legend.title = element_text(size=13),
          axis.title=element_text(size=14),
          legend.direction = "horizontal",
          legend.position = "bottom")
  if(sum(!is.na(New_breaks_labels))==length(New_breaks_labels)){
    pb <- ggplot_build(GG_complete_chi)
    y_breaks <-pb$layout$panel_params[[1]]$y$breaks
    New_breaks_labels<-New_breaks_labels[y_breaks]
    GG_complete_chi<-GG_complete_chi+
      scale_x_continuous(breaks = y_breaks,
                         labels = New_breaks_labels)+
      
      scale_y_continuous(breaks = y_breaks,
                         labels = New_breaks_labels)
  }
  ggsave(filename = Filename,
         plot = GG_complete_chi,
         width=8,height = 6)
  return(GG_complete_chi)
}
#' Evolution with k of the convergence
#' of the angular component for each forcing condition. 
#'
#' @param list_MV_Orig list[str: dataframe]. Multivariate
#' time series where the key corresponds to the variable name.
#' 
#' @param l_name vector[str]. Name of forcing conditions
#' @param vector_k vector[int]. Possible number 
#' of exceedances.
#' @param Ref_RiskF vector[float]. Values of the 
#' compound risk function.
#'
#' @return list.
#' @export
#'
#' @examples
Convgce_Angle_evol<-function(list_MV_Orig,l_name,vector_k,
                             Ref_RiskF){
  
  All_results<-lapply(X = vector_k,
                      Launch_MVconvergence_per_K,
                      l_name=l_name,
                      list_MV_Orig=list_MV_Orig,
                      vect_lg=Ref_RiskF)
  
  return(All_results)
}
#' Value for a given k of the mean absolute 
#' coordinates in a finite-dimensional basis for
#' several forcing conditions. 
#' 
#' @param list_MV_Orig: list[df]. 
#' @param l_name: vect[str].
#' @param k: int. Number of extreme multivariate extreme time series. 
#' @param vect_lg: vect[float]. Value of the risk function (gol). 
#'
#' @return list.
#' @export
#'
#' @examples
Launch_MVconvergence_per_K<-function(list_MV_Orig,l_name,k,
                                     vect_lg){
  
  Order_lg<-order(vect_lg,
                  decreasing = TRUE)
  vect_lg_sort<-sort(vect_lg,
                     decreasing = TRUE)
  Threshold_lg<-vect_lg_sort[k]
  ### Subset of extreme events.
  sub_order<-Order_lg[c(1:k)]
  Df<-list()
  list_conv<-list()
  for(nameV in l_name){
    Sub_Mat<-list_MV_Orig[[nameV]][sub_order,]
    L2_extj<-apply(X = Sub_Mat,MARGIN = 1,
                   FUN = calcul_norm_L2)
    Angle_extj<-t(t(Sub_Mat)%*%diag(L2_extj^(-1)))
    list_conv[[nameV]]<-Function_conv_univ(
      Shape_d =Angle_extj )
  }
  return(list_conv)
}

#' Value for each forcing condition 
#' for a given k of the mean absolute coordinate for several basis functions.
#'
#' @param Shape_d df. Matrix of observations for one forcing condition. 
#'
#' @return Vect[float]. Mean absolute coordinates in 8 basis functions (sin) 
#' for the input df. 
#' @export 
#'
#' @examples
Function_conv_univ<-function(Shape_d){
  
  pas_x<-1/ncol(Shape_d)
  vecteur_temps<-c(1:ncol(Shape_d))/ncol(Shape_d)
  fonction_propre_j<-function(vecteur_temps,j){
    fnct_par_temps<-function(j,t){
      return(sin(2*pi*t*j))
    }
    vecteur_r<-sapply(vecteur_temps,fnct_par_temps,j=j)
    return(vecteur_r)
  }
  list_convergence<-c()
  for(l in c(1:8)){
    fonction_obtenue<-fonction_propre_j(vecteur_temps =vecteur_temps,j=l )
    
    # approx de Rieman --------------------------------------------------------
    coordonnees<-(Shape_d%*%fonction_obtenue)*(pas_x)
    list_convergence<-c(list_convergence,
                         mean(abs(coordonnees)))
  }
  
  return(list_convergence)
}

### base comes from tea::mindist, idea is to better 
### see the variations of the distance 
### to prevent take a not so interesting threshold.
mindist_update<-function (data, ts = 0.15, method = "mad") 
{
  xstat = sort(data, decreasing = TRUE)
  n = length(data)
  T = floor(n * ts)
  i = 1:(n - 1)
  h = (cumsum(log(xstat[i]))/i) - log(xstat[i + 1])
  xstat = sort(data)
  A = matrix(ncol = T - 1, nrow = T - 1)
  for (k in 1:(T - 1)) {
    for (j in 1:(T - 1)) {
      A[k, j] = abs((((k/j) * xstat[n - k + 1]^(1/h[k]))^h[k]) - 
                      xstat[n - j])
    }
  }
  if (method == "mad") {
    M = rowMeans(A)
  }
  if (method == "ks") {
    rowMax <- function(rowData) {
      apply(rowData, MARGIN = c(1), max)
    }
    M = rowMax(A)
  }
  method<-rep(method,length(M))
  df<-cbind.data.frame(1:(T-1),M,method)
  colnames(df)<-c("Nb_k","value_metric","metric")
  return(df)
}
#' Title
#'
#' @param dims_elt_text 
#' @param Vectors_HTAIL 
#' @param Vect_k 
#' @param q 
#' @param root_export 
#' @param YLIM_MV 
#'
#' @return
#' @export
#'
#' @examples
Run_diagnostics_Gamma_G<-function(dims_elt_text,Vectors_HTAIL,
                                  Vect_k,q,root_export,YLIM_MV){
  
  Mat_combs<-t(utils::combn(x = c(1:ncol(Vectors_HTAIL)),
                            m = 2))
  list_gg<-list()
  Names<-colnames(Vectors_HTAIL)
  for(j in c(1:nrow(Mat_combs))){
    pair<-Mat_combs[j,]
    pair_names<-Names[pair]
    Sub_vector<-Vectors_HTAIL[,pair]
    Minj<-apply(X = Sub_vector,MARGIN = 1,FUN = min)
    GG<-Graphics_estimators_gamma(series = Minj,
                                  vect_k = Vect_k,
                                  Title_graphic = " ",
                                  dims_elt_text = dims_elt_text,
                                  SELECT_estim = c("ML_gamma"))
    DF_j<-GG$data
    DF_j$pair_vars<-rep(paste0("Min(",pair_names[1],",",
                               pair_names[2],")"),nrow(DF_j))
    list_gg[[j]]<-DF_j
  }
  Df_all<-do.call(rbind.data.frame,list_gg)
  COLS_chosen<-c("gam_ref"="red",
                 "ML_gamma"="blue",
                 "confidence_band"="darkblue",
                 "threshold"="red")
  LINETYPE_chosen<-c("gam_ref"=3,
                     "Hill_gamma"=2,
                     "ML_gamma"=4,
                     "confidence_band"=5,
                     "threshold"=6)
  Df_all$pair_vars<-sapply(Df_all$pair_vars,
                           FUN =function(x){
                             Fct_correct_name(x = x,target = "Surcote",
                                              replacement = "Surge")})
  GG_shape_AD<-ggplot2::ggplot(data=Df_all,aes(x=number_excesses,y =gamma_estimed,color=source,
                                               group=interaction(source),
                                               linetype=source))+
    geom_line()+
    facet_wrap(~pair_vars)+
    ylab(expression(gamma))+
    xlab("Number of exceedances")+
    ylim(YLIM_MV)+
    geom_ribbon(mapping = aes(ymin=bound_inf,ymax=bound_sup,
                              col="confidence_band",
                              linetype = "confidence_band"),alpha=0.15,
                fill="grey")+
    geom_hline(aes(yintercept=Shape_chosen,col="gam_ref",
                   linetype="gam_ref"))+
    labs(col="Legend", linetype = "Legend")
  LABELS_shape<-names(COLS_chosen)
  LEGEND_shape<-sapply(LABELS_shape,
                       function(x){
                         if(x=="gam_ref"){
                           return(latex2exp::TeX("$\\gamma_{0}$"))
                         }else{
                           return(x)
                         }
                       })
  GG_shape_new<-GG_shape_AD+
    scale_linetype_manual(values = LINETYPE_chosen,
                          labels=LEGEND_shape)+
    scale_color_manual(values = COLS_chosen,
                       labels=LEGEND_shape)+
    theme_bw()+
    theme(axis.title=element_text(size=dims_elt_text[1]),
          legend.title = element_text(size=dims_elt_text[2]),
          legend.text=element_text(size=dims_elt_text[3]),
          axis.text=element_text(size=dims_elt_text[4]))+
    theme(legend.position="bottom",
          legend.direction = "horizontal")
  ggsave(filename = paste0(root_export,"/ADiag_shape_d=",ncol(Vectors_HTAIL),
                           ".png"),
         plot=GG_shape_new,width = 8,height = 6)
  return(NA)
  
  
}

#' Title
#'
#' @param method_corr 
#' @param df1 
#' @param Intersect_times 
#' @param df2 
#' @param mat_corr 
#'
#' @return
#' @export
#'
#' @examples
Fct_correlations<-function(method_corr,df1,
                           Intersect_times,df2=NA,
                           mat_corr=TRUE){
  
  if(is.null(dim(df2))){
    Na_found<-(!is.na(df2))
    Cond<-sum(Na_found)==length(df2)
  }else{
    Na_found<-(is.na(df2)) 
    Test<-sum(colSums(Na_found))
    Cond<-Test==0
  }
  if(Cond){
    vect_r<-sapply(X = Intersect_times,
                   FUN = function(j,x_1,x_2){
                     series_1<-x_1[,j]
                     
                     series_2<-x_2[,j]
                     return(cor(x = series_1,y = series_2,method = method_corr))
                   },x_1=df1,x_2=df2)
    return(vect_r)
  }else{
    d<-ncol(df1)
    NAMES<-colnames(df1)
    if(mat_corr){
      Mat_corr<-cor(df1,method = method_corr)
      return(Mat_corr)
    }else{
      Matrix_cases<-t(utils::combn(x = c(1:d),m = 2))
      list_corr<-list()
      list_pairs<-list()
      for(j in c(1:nrow(Matrix_cases))){
        Pair<-Matrix_cases[j,]
        Pair_used<-NAMES[Pair]
        Sub_df<-df1[,Pair_used]
        value_corr<-cor(x = Sub_df[,1],y = Sub_df[,2],
                        method = method_corr)
        list_corr[[j]]<-value_corr
        list_pairs[[j]]<-paste0("(",Pair_used[1],",",
                                Pair_used[2],")")
      }
      return(data.frame("value"=unlist(list_corr),
                        "pair"=unlist(list_pairs)))
    }
  }
  
}
#' Title
#'
#' @param k_end 
#' @param bandwidth_h 
#' @param series_orig 
#'
#' @return
#' @export
#'
#' @examples
Window_std_per_threshold<-function(k_end,bandwidth_h,series_orig){
  
  Inds_taken<-k_end-bandwidth_h
  ### take results found with higher threshold--> go backwards
  sub_series<-series_orig[Inds_taken:k_end]
  return(sd(sub_series))
}
#' Title
#'
#' @param series_sorted_stats 
#' @param k_max 
#' @param bandwidth_h 
#'
#' @return
#' @export
#'
#' @examples
Stability_criterion_choice_threshold<-function(series_sorted_stats,k_max,
                                               bandwidth_h){
  
  ind_beg<-bandwidth_h
  Inds_candidates<-c(ind_beg:k_max)
  Vect_std_window<-sapply(Inds_candidates,
                          FUN = Window_std_per_threshold,
                          series_orig=series_sorted_stats,
                          bandwidth_h=bandwidth_h)
  ### Detect local minimal for the std vector
  Differences<-diff(Vect_std_window)
  L_end<-length(Differences)-1
  ### The previous index gives a higher value
  Index_1<-which(Differences[1:L_end]<0)+1
  ### The next index gives a higher value
  d2<-Differences
  Index_2<-which(d2>=0)
  ind_min_local<-intersect(Index_1,Index_2)
  
  ### Determine the local minimum giving a value smaller
  ### than the mean
  ref_val<-mean(Vect_std_window)
  All_candidates<-Vect_std_window[ind_min_local]
  sub_ind_candidates<-which(All_candidates<ref_val)
  index_special<-ind_min_local[sub_ind_candidates]
  
  ### Since we go back to the statistic--> we need to
  ### recover the original index--> add bandwidth_h
  Anchor<-index_special[1]+bandwidth_h
  end_anchor<-index_special[1]
  Inds_sub<-c(Anchor:end_anchor)
  ### select the observations close to the chosen beta
  Sub_stat<-series_sorted_stats[Inds_sub]
  
  ### pick among these observations the one closer to the median
  Med_sub_stat<-median(Sub_stat)
  Ind_beta_star<-which.min(abs(Sub_stat-Med_sub_stat))
  beta_star<-Inds_sub[Ind_beta_star]
  ### (graphics) recover the associated standard deviation (if available)
  newbeta_star<-ifelse(beta_star>bandwidth_h,
                       yes = beta_star-bandwidth_h,
                       no = NA)
  print(c(beta_star,newbeta_star))
  STAT_beta<-ifelse(beta_star>bandwidth_h,
                    yes=Vect_std_window[newbeta_star],
                    no=NA)
  return(list("evol_std"=Vect_std_window,
              "beta_found"=newbeta_star,
              "candidates"=index_special,
              "stat_candidates"=Vect_std_window[index_special],
              "stat_beta"=STAT_beta))
}
Concatenate_entries_same_key<-function(list_,key_target){
  vector_results<-c(sapply(which(names(list_)==key_target),
                           FUN=function(x){
                             return(list_[[x]])
                           }))
  return(vector_results)
}

#' Title
#'
#' @param gamma 
#' @param sigma_EGPD 
#' @param kappa 
#' @param order_quantiles 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
Ratio_log_model1<-function(gamma,sigma_EGPD,order_quantiles,data){
  
  Xfound<-as.numeric(quantile(data,order_quantiles))
  logemp_fonction<-log(1-(1+gamma*(Xfound/sigma_EGPD))^(-1/gamma))/log(order_quantiles)
  return(logemp_fonction)
}
RL_EGPD_model1<-function(gamma,sigma_EGPD,order_quantiles,Kappa){
  return(mev::qextgp(p = order_quantiles,kappa =Kappa,
                     sigma = sigma_EGPD,xi = gamma,
                     type = 1))
}
