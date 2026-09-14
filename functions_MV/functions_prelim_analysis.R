
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
  GG_complete_AD_AI<-ggplot(All_AD_AI, aes(x = Col, y = Row, fill = Value)) 
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
#' object.
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
    GGdefault<-ggplot(df, aes(x = Col, y = Row, fill = Value)) +
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
#' Chi_bar_measure_analysis_pair
#'
#' @param Matrix_df_unif dataframe[float]. Univariate
#' matrix with uniform margins
#' @param Order_quantile float. Quantile in the uniform scale.
#' @param t int. Coordinate index to analyse. 
#' @param s int. Other coordinate index to analyse. 
#'
#' @return
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

#' Chi_bar_measure_analysis
#'
#' @param Matrix_df_unif dataframe[float]. Univariate
#' matrix with uniform margins
#' @param Order_quantile float. Quantile in the uniform scale.
#' @param Filename string. Name of the exported graphic. 
#' @param Name_main String (NA by default). Name of the variable 
#' analysed (used for the title of ggplot2 object). 
#'
#' @return Ggplot object summarising the asymptotic 
#' dependencies between the components of the input matrix.
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
    #return(All)
    return(t(Mat_Xibar_found))
  }else{
    mat <-All
    rownames(mat) <- as.character(c(1:ncol(mat)))
    colnames(mat)<- as.character(c(1:ncol(mat)))
    
    # Convert matrix to long format
    df <- melt(mat)
    require(viridis)
    input_for_tex<-paste0("$\\bar{\\chi}$ matrix of ",Name_main)
    Name_main_modif<-latex2exp::TeX(input_for_tex)
    colnames(df) <- c("Row", "Col", "Value")
    GGdefault<-ggplot(df, aes(x = Col, y = Row, fill = Value)) +
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
  GG_complete_chi<-ggplot(All_chi, aes(x = Col, y = Row, fill = Value)) +
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
#' @param liste_MV_Orig list[str: dataframe]. Multivariate
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
Convgce_Angle_evol<-function(liste_MV_Orig,l_name,vector_k,
                             Ref_RiskF){
  
  All_results<-lapply(X = vector_k,
                      Launch_MVconvergence_per_K,
                      l_name=l_name,
                      liste_MV_Orig=liste_MV_Orig,
                      vect_lg=Ref_RiskF)
  
  return(All_results)
}
#' Value for a given k of the mean absolute 
#' coordinates in a finite-dimensional basis for
#' several forcing conditions. 
#' 
#' @param liste_MV_Orig: list[df]. 
#' @param l_name: vect[str].
#' @param k: int. Number of extreme multivariate extreme time series. 
#' @param vect_lg: vect[float]. Value of the risk function (gol). 
#'
#' @return list.
#' @export
#'
#' @examples
Launch_MVconvergence_per_K<-function(liste_MV_Orig,l_name,k,
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
    Sub_Mat<-liste_MV_Orig[[nameV]][sub_order,]
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
  LISTE_convergence<-c()
  for(l in c(1:8)){
    fonction_obtenue<-fonction_propre_j(vecteur_temps =vecteur_temps,j=l )
    
    # approx de Rieman --------------------------------------------------------
    coordonnees<-(Shape_d%*%fonction_obtenue)*(pas_x)
    LISTE_convergence<-c(LISTE_convergence,
                         mean(abs(coordonnees)))
  }
  
  return(LISTE_convergence)
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