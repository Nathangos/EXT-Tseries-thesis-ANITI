rm(list=ls())

##################
set.seed(133)
source("../fonctions/fonctions_perfs_ML.R")
source("../fonctions/fonctions.R")
source("functions_MV/Fonctions_consistency.R")
source("functions_MV/function_estim_simul_gPareto.R")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(extRemes)
require(ggside)
require(patchwork)
require(viridis)
# Import residuals --------------------------------------------------------
cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")
variable<-"Surcote"

# choice of the risk function ---------------------------------------------
Name_riskF<-"Weighted_NQ"
Name_riskF2<-Name_riskF
#Name_riskF<-"max"
q<-(20)
W_chosen<-c(0.3,0.7)
if(Name_riskF=="sum"){
  Risk_f<-function(x){return(sum(x))}
  wFun<-weightFun_sum
  dwFun<-DweightFun_sum
  
}
if(Name_riskF=="max"){
  Risk_f<-function(x){return(max(x))}
}
if(Name_riskF=="sum_penalized"){
  Risk_f<-function(x){return(sum_penalized(x = x,vect_w = c(0.7,0.3)))}
  wFun<-function(x,u){
    weightFun_sum_penalized(x = x,u = u,
                            vect_w = W_chosen)
  }
  dwFun<-function(x,u){
    DweightFun_sum_penalized(x = x,u = u,
                             vect_w = W_chosen)
  }
}

if(Name_riskF=="Weighted_NQ"){
  Name_riskF2<-paste0(Name_riskF,"_q=",q)
  Risk_f<-function(x){
    return(Norm_q(x=x,q=q,
                  weights_nq = c(0.7,0.3)))
  }
  wFun<-function(x,u){
    return(weightFun_nq(x = x,u = u,
                        q = q,
                        weights_nq =c(0.7,0.3)))
  }
  dwFun<-function(x,u){
    return(DweightFun_nq(x = x,u = u,
                         q = q,weights_nq = c(0.7,0.3)))
  }
}
ANALYSIS_for_resid<-TRUE
type_donnees<-"HIVER"
L_nameV<-c("Surcote","U")
if(ANALYSIS_for_resid){
  l_root<-"residuals_MV/"
  path <- paste0(getwd(),"/residuals_MV/")
  ROOT_import<-paste0(path,Name_riskF2,"/")
  CPLMT<-"_resid"
  path_graph<-paste0(getwd(),"/graphiques_MV/resid")
 
}else{
  l_root<-"../ss_tend/"
  if(type_donnees=="HIVER"){
    l_root<-paste0(l_root,"HIVER/")
  }
  path <- paste0(getwd(),"/Generated_extremes/")
  ROOT_import<-paste0(path,Name_riskF2,"/")
  CPLMT<-"_data"
  path_graph<-paste0(getwd(),"/graphiques_MV/data")
}
path_graph
## Create the new sub directory if needed
try({
  dir.create(file.path(path_graph))
})
## Create for the specific risk function used the sub directory used.
newfolder <- Name_riskF2
try({
  dir.create(file.path(path_graph, newfolder))
})
Updated<-paste0(path_graph,"/",Name_riskF2)
ROOT_export_graph<-paste0(path_graph,"/",
                          Name_riskF2,"/")
ROOT_export_graph
try({
  TEND_path<-paste0(ROOT_export_graph,"Tend")
  dir.create(file.path(TEND_path))
  Theta_path<-paste0(ROOT_export_graph,"Theta")
  dir.create(file.path(Theta_path))
  Lg_path<-paste0(ROOT_export_graph,"RisKfunctions")
  dir.create(file.path(Lg_path))
  PCA_path<-paste0(ROOT_export_graph,"PCA")
  dir.create(file.path(PCA_path))
  Extremo_path<-paste0(ROOT_export_graph,"extremes")
  dir.create(file.path(Extremo_path))
  for(name_V in L_nameV){
    newfolder<-name_V
    try({
      dir.create(file.path(TEND_path, newfolder))
      dir.create(file.path(Theta_path, newfolder))
      dir.create(file.path(PCA_path, newfolder))
      dir.create(file.path(Extremo_path, newfolder))
    })
  }
})
if(ANALYSIS_for_resid){
  Simuls<-read.csv(file=paste0(ROOT_import,variable,"_simul.csv"))[,2:38]
  Reals<-read.csv(file=paste0(ROOT_import,variable,"_obs_ext.csv"))[,2:38]
}else{
  Simuls<-read.csv(file=paste0(ROOT_import,variable,"_fsim.csv"))[,2:38]
  Reals<-read.csv(file=paste0(ROOT_import,variable,"_fdata.csv"))[,2:38]
}

colnames(Simuls)<-c(1:ncol(Simuls))
colnames(Reals)<-colnames(Simuls)
Nsim<-nrow(Simuls)


# PCA coordinates ---------------------------------------------------------
if(ANALYSIS_for_resid){
  N_sum<-0
  Rank_Cop<-1
  LIST_GG<-list()
  Fich_1<-read.csv(file = paste0(ROOT_import,"coord_data.csv"))
  Fich_2<-read.csv(file = paste0(ROOT_import,"coord_simul.csv"))
  Fich_1<-Fich_1[,c(2:ncol(Fich_1))]
  Fich_2<-Fich_2[,c(2:ncol(Fich_2))]
  colnames(Fich_1)<-c(1:ncol(Fich_1))
  colnames(Fich_2)<-colnames(Fich_1)
  Combinaisons<-t(utils::combn(x = c(1:ncol(Fich_1)),m = 2))
  for(j in c(1:nrow(Combinaisons))){
    index_i<-Combinaisons[j,1]
    index_j<-Combinaisons[j,2]
    Pair_<-c(index_i,index_j)
    df_combo<-rbind.data.frame(Fich_1[,Pair_],
                               Fich_2[,Pair_])
    colnames(df_combo)<-c("V1","V2")
    df_combo$Legend<-c(rep("data",nrow(Fich_1)),
                       rep("simulations",nrow(Fich_2)))
    cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")
    Lab_i<-expression("C"[index_i])
    Lab_i_eval<-do.call("substitute", list( Lab_i[[1]], 
            list(index_i =index_i)))
    Lab_j<-expression("C"[index_j])
    Lab_j_eval<-do.call("substitute", list( Lab_j[[1]], 
                                            list(index_j =index_j)))
    GG1<-ggplot(data=df_combo,aes(x=V1,y=V2,colour=Legend,shape=Legend),
    )+
      geom_point(aes(size=Legend))+
      geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
      geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
      xlab(Lab_i_eval)+ylab(Lab_j_eval)+scale_color_manual(values=cols_)+
      scale_fill_manual(values=cols_)+
      scale_shape_manual(values = c("simulations"=17,"data"=19))+
      scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
      guides(fill="none")+
      theme(axis.title=element_text(size=15),
            legend.text=element_text(size=10))
    LIST_GG[[j]]<-GG1
    N_sum<-N_sum+1
    if(N_sum==3){
      #blank <- ggplot(width=0.5) + theme_void()
      obj1<-LIST_GG[[3*(Rank_Cop-1)+1]] 
      obj2<-LIST_GG[[3*(Rank_Cop-1)+2]]
      obj3<-LIST_GG[[3*(Rank_Cop)]]
      top_row<-(obj1+obj2)
      Second_row<-plot_spacer() + obj3 + plot_spacer()+
        plot_layout(widths = c(0.5, 1,0.5))
      All_PCA_coords<-top_row/Second_row
      Name_graph<-paste0("/Copula_model_vs_coords_",Rank_Cop,".png")
      ggsave(filename =paste0(Theta_path,Name_graph),
             plot =All_PCA_coords,width = 12,height = 9)
      Rank_Cop<-Rank_Cop+1
      N_sum<-0
    }
  }
}
# Plots with more than 1 forcing condition --------------------------------
##########
Link_RF<-paste0("residuals_MV/",Name_riskF2,"/")
Unit_variable<-list("Surcote"="(m)","U"="(m/s)")
list_simul<-list()
list_obs<-list()
Obs_taken<-sample(c(1:nrow(Reals)),size = 100,
                  replace = FALSE)
Simul_taken<-sample(c(1:Nsim),size = 100,
                    replace = FALSE)
Simul_taken
# Params for extremo computation 
# Extremogram -------------------------------------------------------------
Matrice_couples<-list()
L<-ncol(Reals)
z<-1
vecteur_distances<-c()
for(j in 1:L){
  #condition imposée sur le deuxième temps.
  valeurs_t_plus_h<-j:L
  for (i in valeurs_t_plus_h){
    Matrice_couples[[z]]<-c(i,j)
    vecteur_distances<-c(vecteur_distances,abs(i-j))
    z<-z+1
  }
}
q_chosen<-0.70
Times_chosen<-c(19,13,25,31)
periods_years<-c(2,5,10,20,50,80,
                 100,120,200,400)
NB_years<-37
Qtile_per_year<-0.90
for(NV in L_nameV){
  Z<-which(L_nameV==NV)
  if(ANALYSIS_for_resid){
    Simuls_1<-read.csv(file=paste0(ROOT_import,NV,"_simul.csv"))[,2:38]
    Reals_1<-read.csv(file=paste0(ROOT_import,NV,"_obs_ext.csv"))[,2:38]
  }else{
    Simuls_1<-read.csv(file=paste0(ROOT_import,NV,"_fsim.csv"))[,2:38]
    Reals_1<-read.csv(file=paste0(ROOT_import,NV,"_fdata.csv"))[,2:38]
  }
  colnames(Simuls_1)<-c(1:ncol(Simuls_1))
  colnames(Reals_1)<-colnames(Simuls_1)
  list_simul[[NV]]<-Simuls_1
  list_obs[[NV]]<-Reals_1

  # Violin plots/Boxplots ------------------------------------------------------------
  ##########
  Melt_sim_violin<-melt(t(Simuls_1))
  Melt_data_violin<-melt(t(Reals_1))
  colnames(Melt_sim_violin)<-c("time","id","value")
  colnames(Melt_data_violin)<-c("time","id","value")
  Melt_whole<-rbind.data.frame(Melt_data_violin,
                               Melt_sim_violin)
  Melt_whole$origin<-c(rep("data",nrow(Melt_data_violin)),
                       rep("simulations",nrow(Melt_sim_violin)))
  Melt_whole$time<-as.character(( Melt_whole$time-19)/6)
  GG_violin<-ggplot(data=Melt_whole,aes(x = time,y=value))+
    geom_boxplot(width=0.8,aes(fill=origin))+
    facet_wrap(~origin)+
    labs(fill="Legend")+
    ylab(paste0("value ",Unit_variable[[NV]]))+
    theme(axis.text.x = element_blank())+
    scale_fill_manual(values=cols_)
  ggsave(filename = paste0(TEND_path,"/",NV,"/",NV,"_boxplot_sim_obs",CPLMT,".png"),
         plot = GG_violin,
         width = 8,height = 6)
  # QQ evol -----------------------------------------------------------------
  #####
  # Simul -------------------------------------------------------------------
  estimator_95_simul<-apply(Simuls_1,MARGIN =2,
                            function(x){return(as.numeric(quantile(x,0.95)))})
  estimator_05_simul<-apply(Simuls_1,MARGIN =2,
                            function(x){return(as.numeric(quantile(x,0.05)))})
  estimator_50_simul<-apply(Simuls_1,MARGIN =2,
                            function(x){return(as.numeric(quantile(x,0.50)))})
  estimator_975_simul<-apply(Simuls_1,MARGIN =2,
                             function(x){return(as.numeric(quantile(x,0.975)))})
  
  ### data estimator
  estimator_05<-apply(Reals_1,MARGIN =2,
                      function(x){return(as.numeric(quantile(x,0.05)))})
  estimator_95<-apply(Reals_1,MARGIN =2,
                      function(x){return(as.numeric(quantile(x,0.95)))})
  estimator_975<-apply(Reals_1,MARGIN =2,
                       function(x){return(as.numeric(quantile(x,0.975)))})
  estimator_05<-apply(Reals_1,MARGIN =2,
                      function(x){return(as.numeric(quantile(x,0.05)))})
  estimator_50<-apply(Reals_1,MARGIN =2,
                      function(x){return(as.numeric(quantile(x,0.50)))})
  
  # Obs-confidence -----------------------------------------------------------
  N_rep<-500
  B<-nrow(Reals_1)
  N_rep<-500
  lq<-c(0.05,0.50,0.95,0.975)
  Tend<-replicate(N_rep,Resamples_tendencies(B = B,extremes_indus= Reals_1,
                        list_Q =lq))
  
  bound_minus05<-apply(Tend[1,,],MARGIN = 1,
                       FUN = function(x){return(quantile(x,0.025))})
  bound_plus05<-apply(Tend[1,,],MARGIN = 1,
                      FUN = function(x){return(quantile(x,0.975))})
  
  ##########
  bound_minus50<-apply(Tend[2,,],MARGIN = 1,
                       FUN = function(x){return(quantile(x,0.025))})
  bound_plus50<-apply(Tend[2,,],MARGIN = 1,
                      FUN = function(x){return(quantile(x,0.975))})
  
  ########
  bound_minus95<-apply(Tend[3,,],MARGIN = 1,
                       FUN = function(x){return(quantile(x,0.025))})
  bound_plus95<-apply(Tend[3,,],MARGIN = 1,
                      FUN = function(x){return(quantile(x,0.975))})
  
  ######
  bound_minus975<-apply(Tend[4,,],MARGIN = 1,
                        FUN = function(x){return(quantile(x,0.025))})
  bound_plus975<-apply(Tend[4,,],MARGIN = 1,
                       FUN = function(x){return(quantile(x,0.975))})
  
  # Graphics ----------------------------------------------------------------
  df<-cbind.data.frame(c(estimator_05,estimator_50,estimator_95,estimator_975),
                       c(bound_minus05,bound_minus50,
                         bound_minus95,bound_minus975),
                       c(bound_plus05,bound_plus50,
                         bound_plus95,bound_plus975))
  colnames(df)<-c("estimator","bound_minus","bound_plus")
  
  df$estimator_simul<-c(estimator_05_simul,estimator_50_simul
                        ,estimator_95_simul,estimator_975_simul)
  df$Time<-rep(c(1:ncol(Reals_1)),4)
  df$percent<-c(rep("Q05",ncol(Reals_1)),
                rep("Q50",ncol(Reals_1)),
                rep("Q95",ncol(Reals_1)), 
                rep("Q975",ncol(Reals_1)))
  GG_percent<-ggplot(data=df,aes(x=Time,y=estimator,group=interaction(percent),col="data"))+
    facet_wrap(~percent,scales="free_y")+
    geom_line()+
    geom_point()+
    geom_ribbon(mapping = aes(ymin=bound_minus,ymax=bound_plus,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    geom_line(aes(x=Time,y=estimator_simul,col="simulations"),linetype=2)+
    theme(axis.title=element_text(size=15))+
    ylab(paste0("value ",Unit_variable[[NV]]))+
    scale_color_manual(values=cols_)+
    xlab("Time")+
    labs(col="Legend")
  ggsave(filename = paste0(TEND_path,"/",NV,"/QQ_evol",NV,CPLMT,".png"),
         plot = GG_percent,width = 8,height = 6)
  # Plots  ------------------------------------------------------------------
  ###########
  Q975<-apply(X = Reals_1,
              MARGIN =2,FUN = function(x){
                return(quantile(x,0.975))
              })
  Q025<-apply(X = Reals_1,
              MARGIN = 2,FUN=function(x){
                return(quantile(x,0.025))
              })
  Q975_sim<-apply(X = Simuls_1,
              MARGIN = 2,FUN=function(x){
                return(quantile(x,0.975))
              })
  Q025_sim<-apply(X = Simuls_1,
              MARGIN =2,FUN= function(x){
                return(quantile(x,0.025))
              })
  Means_var<-colMeans(Reals_1)
  Means_var_sim<-colMeans(Simuls_1)
  Obs_shown<-Reals_1[Obs_taken,]
  Simul_shown<-Simuls_1[Simul_taken,]
  Melt_data<-melt(t(Obs_shown))
  colnames(Melt_data)<-c("time","id","value")
  Melt_data$mean<-Means_var[Melt_data$time]
  Melt_data$Qinf<-Q025[Melt_data$time]
  Melt_data$Qsup<-Q975[Melt_data$time]
  
  Melt_sim<-melt(t(Simul_shown))
  colnames(Melt_sim)<-c("time","id","value")
  Melt_sim$mean<-Means_var_sim[Melt_sim$time]
  Melt_sim$Qinf<-Q025_sim[Melt_sim$time]
  Melt_sim$Qsup<-Q975_sim[Melt_sim$time]
  Melt_all<-rbind.data.frame(
    Melt_data,Melt_sim
  )
  Melt_all$origin<-c(rep("data",nrow(Melt_data)),
                     rep("simulations",nrow(Melt_sim)))
  Melt_all$time<-(Melt_all$time-19)/6

  min_y<-min(min(Melt_all$Qinf),min(Melt_all$value))
  max_y<-max(max(Melt_all$Qsup),max(Melt_all$value))
  GG_simul_vs_obs<-ggplot(data=Melt_all,aes(x=time,y=value
                           ,group=interaction(id),
                           col=id))+
    guides(col="none")+
    scale_color_viridis() +
    geom_line()+
    facet_wrap(~origin)+
    geom_ribbon(aes(ymin=Qinf,ymax=Qsup,linetype="confidence_band"),
                alpha=0.02,fill="grey",col="darkblue")+
    geom_line(aes(x = time,y=mean,linetype="mean"),
              col="red",size=1.1)+
    geom_line(alpha=0.7)+
    scale_linetype_manual("Legend",values=c("confidence_band"=2,
                                            "mean"=5))+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    xlab(" ")+
    ylab(paste0("value ",Unit_variable[[NV]]))+
    ylim(c(min_y,max_y))+
    theme(axis.title=element_text(size=25),
          legend.text=element_text(size=14),
          legend.title = element_text(size=15),
          axis.text = element_text(size=12),
          strip.text = element_text(size = 12))
  ggsave(filename = paste0(TEND_path,"/",NV,"/Tend_",NV,"_simul_obs",CPLMT,".png"),
         plot = GG_simul_vs_obs,
         width = 8,height = 6)
  
  # Extremogram 
  ########
  Tau1<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_chosen,
    x=Simuls_1)
  Tau2<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_chosen,
    x=Reals_1)
  Extremogram_simulations<-empirical_extremogram(Matrix_couples =Matrice_couples,Tau = Tau1,
                                                 inds_select= Simuls_1)
  Extremogram_data<-empirical_extremogram(Matrix_couples=Matrice_couples,
                                          inds_select = Reals_1,Tau = Tau2)
  N_rep<-500
  B<-nrow(Reals_1)
  Result<-replicate(N_rep,fnct_estim_extremo_resample(B = B,inds_ext = Reals_1,
                                                       Tau = Tau2,Matrix_couples =  Matrice_couples,
                                                       vector_distances=vecteur_distances))
  fnct_quantile<-function(x,q){
    M<-apply(X = x,MARGIN = 2,FUN = function(X){return(as.numeric(quantile(X,q,
                                                      )))})
    return(M)
  }
  Q05<-fnct_quantile(x = t(Result),q=0.025)
  Q95<-fnct_quantile(x = t(Result),q=0.975)
  df<-cbind.data.frame(Extremogram_simulations,Extremogram_data,vecteur_distances)
  colnames(df)<-c("simulation_value","data_value","delta")
  result_delta<-df %>% group_by(delta) %>% summarise(val_data=mean(data_value),
                                                     val_simul=mean(simulation_value))
  result_delta$bound_inf<-Q05
  result_delta$bound_sup<-Q95
  extremo_<-as.data.frame(result_delta[,c("val_data","val_simul","bound_inf","bound_sup")])
  extremo_$Time<-c(1:nrow(extremo_))
  GG_extremo<-ggplot(extremo_,aes(x=Time,y=val_data,col="data"))+
    geom_line()+
    geom_point()+
    geom_line(aes(y=val_simul,col="simulations"))+
    geom_point(aes(y=val_simul,col="simulations"),pch=2)+
    geom_ribbon(mapping = aes(ymin=bound_inf,ymax=bound_sup,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    scale_color_manual(values=cols_)+
    theme(axis.title=element_text(size=15))+
    labs(col="Legend")+
    xlab("Lag h")+
    ylab("value")
  ggsave(filename = paste0(Extremo_path,"/",NV,"/",NV,"_extremo_univ",CPLMT,".png"),
         plot = GG_extremo,
         width = 8,height = 6)
  ### Return levels
  LIST_GG_RL<-list()
  for(ind_t in c(1:length(Times_chosen))){
    t<-Times_chosen[ind_t]
    if(ANALYSIS_for_resid){
      link_data<-paste0(l_root,NV,"_residuals.csv")
      Data_obs<-read.csv(file=link_data)[,2:38]
      colnames(Data_obs)<-c(1:37)
      rownames(Data_obs)<-c(1:nrow(Data_obs))
    }else{
      link_data<-paste0(l_root,NV,"_ss_tend.csv")
      Data_obs<-read.csv(file=link_data)[,2:38]
      colnames(Data_obs)<-c(1:37)
    }
    Indexes_exceed<-read.csv(file=paste0(Link_RF,
                                         "ext_series_chosen_risks.csv"))
    Indexes_exceed<-Indexes_exceed[,3]
    series_simul<-Simuls_1[,t]
    series_t<-Data_obs[,t]
    threshold<-quantile(series_t,Qtile_per_year)
    npy<-length(series_t)/NB_years
    series_ext<-subset(series_t,series_t>threshold)
    P_lim<-3*NB_years
    periods_years<-periods_years[which(periods_years<P_lim)]
    P<-1-(npy*periods_years)^(-1)
    Y_graph_lim<-c(0.80*threshold,max(max(series_simul),
                                      max(series_t))*1.20)
    GG_rel<-RL_ggplot_cond_ext(series = series_t,seuil = threshold,
                               period_years = periods_years,
                               NPY = npy,plus_simul = TRUE,
                               series_simul = series_simul, 
                               titre = paste0("Return level (ML) ",comment, " the tidal peak for ",Nom_graph), 
                               nom_variable =NV,cols_ggplot = cols_, 
                               ylim_opt = Y_graph_lim,alpha=0.05,
                               Individus_exts=Indexes_exceed)$GG_plot
    LIST_GG_RL[[ind_t]]<-GG_rel
  }
  NB_divisions<-length(Times_chosen)%/%2
  for(Z in c(1:NB_divisions)){
    Ind1<-(Z-1)*2+1
    Ind2<-Z*2
    Time1<-Times_chosen[Ind1]
    Time2<-Times_chosen[Ind2]
    Name_file<-paste0("rlevel_simul_vs_obs_",NV,"_",
                      Time1,"_",Time2,
                      ".png")
    yleft <- textGrob(paste0("return level (",Unit_variable[[NV]],")"),
                      rot=90,
                      gp = gpar(col = "black", fontsize = 25))
    Xbottom<-textGrob("Period P (years)",
                      gp = gpar(col = "black", fontsize = 25))
    Object_RL_twotimes<-gridExtra::grid.arrange(LIST_GG_RL[[Ind1]],
                        LIST_GG_RL[[Ind2]],ncol=2,
                        left=yleft,bottom=Xbottom)
    ggsave(filename = paste0(Extremo_path,"/",
                             NV,"/",
                             Name_file),
           plot = Object_RL_twotimes,width = 14,
           height = 7)
  }
  
  Z<-Z+1  
  # PCA coordinates ---------------------------------------------------------
  Nb_scores<-3
  L2_simul<-apply(X = Simuls_1,MARGIN = 1,
                  FUN = calcul_norme_L2)
  L2_reality<-apply(X = Reals_1,MARGIN = 1,
                    FUN = calcul_norme_L2)
  POS_sim<-which(L2_simul>0)
  Angle_simul<-t(t(Simuls_1[POS_sim,])%*%diag(L2_simul[POS_sim]^(-1)))
  POS<-which(L2_reality>0)
  Angle_inds<-t(t(Reals_1[POS,])%*%diag(L2_reality[POS]^(-1)))
  PCA_<-FactoMineR::PCA(X = Angle_inds,scale.unit = TRUE,
                        graph = FALSE)
  
  percent_variance<-PCA_$eig[c(1:Nb_scores),2]
  # Calculate mean and standard deviation -----------------------------------
  sd_<-apply(X =Angle_inds,MARGIN = 2,FUN = sd)
  mu_<-colMeans(Angle_inds)
  Angle_simul_std<-scale(Angle_simul,center = mu_,
                         scale = sd_)
  
  # Obtain the coordinates for each dimension -------------------------------
  coords_data<-PCA_$ind$coord[,c(1:Nb_scores)]
  V<-PCA_$svd$V[,c(1:Nb_scores)]
  coords_simul<-Angle_simul_std%*%V
  KS_PCA_coords_<-rep(NA,Nb_scores)
  for(g in c(1:Nb_scores)){
    KS_PCA_coords_[g]<-ks.test(coords_simul[,g],
                               coords_data[,g])$p.value
  }
  df_PCA_simul<-as.data.frame(coords_simul)
  summary(df_PCA_simul)
  colnames(df_PCA_simul)<-sapply(c(1:ncol(df_PCA_simul)), 
                                 function(x){return(paste0("Score_",x))})
  
  df_PCA_data<-as.data.frame(coords_data)
  colnames(df_PCA_data)<-sapply(c(1:ncol(df_PCA_data)), 
                                function(x){return(paste0("Score_",x))})
  df_combo<-rbind.data.frame(df_PCA_data,df_PCA_simul)
  df_combo$type<-c(rep("data",nrow(df_PCA_data)),
                   rep("simulations",nrow(df_PCA_simul)))
  GG0<-ggplot(data=df_PCA_data,aes(x=Score_1,y=Score_2,col="data"))+
    geom_point()+
    geom_point(data=df_PCA_simul,aes(x=Score_1,y=Score_2,col="simulations"),
               pch=17,size=0.75)+
    geom_xsidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
    xlab(paste0("First dimension (",round(percent_variance[1],1),"% of the variance)"))+
    ylab(paste0("Second dimension (",round(percent_variance[2],1),"% of the variance)"))+
    labs(col="Legend")+
    scale_color_manual(values=cols_)+
    scale_fill_manual(values=cols_)+
    guides(fill="none")+
    theme(axis.title=element_text(size=15))
  ggsave(filename = paste0(PCA_path,"/",NV,"/PCA_coords_",NV,CPLMT,".png"),
         width=6,height = 5,
         plot = GG0)
  
}
# Correlation values ------------------------------------------------------
###########
NL2_simul<-sapply(X = list_simul,
       function(x){return(apply(x,MARGIN=1,
                                    FUN = calcul_norme_L2))})
NL2_obs<-sapply(X = list_obs,
                  function(x){return(apply(x,MARGIN=1,
                                           FUN = calcul_norme_L2))})
NL2_obs<-as.data.frame(NL2_obs)
NL2_simul<-as.data.frame(NL2_simul)
colnames(NL2_obs)<-L_nameV
colnames(NL2_simul)<-colnames(NL2_obs)
NL2_obs<-as.data.frame(NL2_obs)
Risk_f_obs1<-apply(X = NL2_obs,MARGIN = 1,FUN = Risk_f)
Risk_f_obs<-cbind.data.frame(Risk_f_obs1,rep("data",
                                             length(Risk_f_obs1)))

NL2_simul<-as.data.frame(NL2_simul)
Sub_simul<-colnames(NL2_obs)

Risk_f_simul1<-apply(X = NL2_simul,MARGIN = 1,FUN = Risk_f)

# QQ plot analysis --------------------------------------------------------
Obj_qq<-qqplot(x = Risk_f_obs1,y = Risk_f_simul1,make.plot = FALSE)
QQdata<-Obj_qq$qdata
GGQQ_RISKFsim_data<-ggplot(data=QQdata,aes(x=x,y=x,col="data"))+
  geom_line()+
  geom_point(aes(x=x,y=y,col="simulations"),
             pch=17)+
  scale_color_manual(values = cols_)+
  labs(col="Legend")+
  theme(axis.title=element_text(size=15))+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")
ggsave(plot = GGQQ_RISKFsim_data,
       filename=paste0(Lg_path,"/QQ_RiskF_data_sim",CPLMT,".png"),
       width=8,height=6)

# Violin plots for each dataset (simulations and observations) ------------------------------------------
#####
Orig_s<-rep("simulations",length(Risk_f_simul1))
Risk_f_simul<-cbind.data.frame(Risk_f_simul1,Orig_s)
colnames(Risk_f_simul)<-c("Risk_f","orig")
colnames(Risk_f_obs)<-colnames(Risk_f_simul)
Risk_f_obs_simul<-rbind.data.frame(Risk_f_obs,Risk_f_simul)
colnames(Risk_f_obs_simul)<-c("value_Risk_f","origin")
GG_violin_Risk_f<-ggplot(data=Risk_f_obs_simul,aes(x=origin,y=value_Risk_f,
                              fill=origin))+
  geom_violin()+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values=cols_)+
  labs(fill="Legend")+ylab("value (m)")
ggsave(plot=GG_violin_Risk_f,
       filename=paste0(Lg_path,"/Violin_RiskF_data_sim",CPLMT,".png"),
       width=8,height=8)

NL2_obs$origin<-rep("data",nrow(NL2_obs))
NL2_simul$origin<-rep("simulations",
                      nrow(NL2_simul))
NL2_all<-rbind.data.frame(NL2_obs,NL2_simul)
NL2_all<-melt(NL2_all)
colnames(NL2_all)<-c("origin","variable",
                     "value")
GG_violin<-ggplot(data=NL2_all,aes(x=origin,y=value,fill=origin))+
  facet_wrap(~variable,scales = "free")+
  geom_violin(alpha=0.6)+
  scale_fill_manual(values=cols_)+
  geom_boxplot(width=0.3)+
  labs(fill="Legend")+
  xlab("origin")+ylab("value (m)")

ggsave(plot=GG_violin,
       filename =paste0(Lg_path,"/Violin_RiskF_data_sim_varwise",CPLMT,".png"),
       width=8,height=8)

N_rep<-500
Mcorr<-"kendall"
B<-nrow(Reals_1)
Conf_int_corr<-replicate(n = N_rep,Resamples_correlations(l_extremes_indus = list_obs,
                       B = B,method_corr = Mcorr,
                       l_names = L_nameV))

# Compute quantile at each time -------------------------------------------
Q_minus<-apply(X = Conf_int_corr,MARGIN = 1,
               function(x){return(quantile(x,0.025))})
Q_plus<-apply(X = Conf_int_corr,MARGIN = 1,
               function(x){return(quantile(x,0.975))})
Estimator<-Fct_correlations(method_corr = Mcorr ,
                 df1 =list_obs[[L_nameV[1]]],
                 df2=list_obs[[L_nameV[2]]])
Estimator_corr_simul<-Fct_correlations(method_corr = Mcorr ,
                                       df1 =list_simul[[L_nameV[1]]],
                                       df2=list_simul[[L_nameV[2]]])
DF<-cbind.data.frame(Q_minus,Q_plus,Estimator ,
                     Estimator_corr_simul)
colnames(DF)<-c("bound_min","bound_plus",
                "estimator","estimator_simul")
DF$Time<-c(1:nrow(DF))
GG_corr<-ggplot(data=DF,aes(x=Time,y=estimator,col="data"))+
  geom_line()+
  geom_ribbon(aes(ymin=bound_min,
                  ymax=bound_plus,col="confidence_band"),alpha=0.15,
  fill="grey", linetype = "dashed")+
  geom_line(aes(y=estimator_simul,col="simulations"),linetype=2)+
  theme(axis.title=element_text(size=15))+
  ylab("correlation coefficient")+
  scale_color_manual(values=cols_)+
  xlab("Time")+
  labs(col="Legend")
ggsave(filename = paste0(ROOT_export_graph,"Coor_",Mcorr,CPLMT,".png"),
       plot = GG_corr,
       width=8,height=8)

# Bivariate relation ------------------------------------------------------
L_simul<-NULL
L_obs<-NULL
L_L2_simul<-list()
L_L2_obs<-list()
for(name_variable_analyze in L_nameV){
  if(ANALYSIS_for_resid){
    Simuls_v<-read.csv(file=paste0(ROOT_import,name_variable_analyze,
                                   "_simul.csv"))[,2:38]
    Reals_v<-read.csv(file=paste0(ROOT_import,name_variable_analyze
                                  ,"_obs_ext.csv"))[,2:38]
  }else{
    Simuls_v<-read.csv(file=paste0(ROOT_import,name_variable_analyze,
                                   "_fsim.csv"))[,2:38]
    Reals_v<-read.csv(file=paste0(ROOT_import,name_variable_analyze
                                  ,"_fdata.csv"))[,2:38]
  }
  
  
  colnames(Simuls_v)<-c(1:ncol(Simuls_v))
  colnames(Reals_v)<-colnames(Simuls_v)
  L_L2_simul[[name_variable_analyze]]<-apply(X = Simuls_v,
          MARGIN = 1,FUN = calcul_norme_L2)
  L_L2_obs[[name_variable_analyze]]<-apply(X = Reals_v,
          MARGIN = 1,FUN = calcul_norme_L2)
  L_simul[[name_variable_analyze]]<-Simuls_v
  L_obs[[name_variable_analyze]]<-Reals_v
}
DF_L2_simul<-as.data.frame(L_L2_simul)
DF_L2_obs<-as.data.frame(L_L2_obs)
Whole_L2<-rbind.data.frame(DF_L2_simul,DF_L2_obs)
Whole_L2$Legend<-c(rep("simulations",nrow(DF_L2_simul)),
                   rep("data",nrow(DF_L2_obs)))
if(length(L_nameV)>2){
  Combn_rf<-utils::combn(c(1:length(L_nameV)),m = 2)
}else{
  Combn_rf<-matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)
}
apply(X = Combn_rf,MARGIN = 1,
      function(x){
        name_i<-L_nameV[x[1]]
        name_j<-L_nameV[x[2]]
        Sub_whole_j<-Whole_L2[,c(name_i,name_j,"Legend")]
        colnames(Sub_whole_j)<-c("V1","V2","Legend")
        GG_RF_j<-ggplot(data=Sub_whole_j,
               mapping = aes(x=V1,y=V2,colour=Legend,shape=Legend))+
          geom_point(aes(size=Legend))+
          geom_xsidedensity(data=Sub_whole_j,aes(fill=Legend), alpha = 0.5)+
          geom_ysidedensity(data=Sub_whole_j,aes(fill=Legend), alpha = 0.5)+
          scale_color_manual(values=cols_)+
          scale_fill_manual(values=cols_)+
          scale_shape_manual(values = c("simulations"=17,"data"=19))+
          scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
          guides(fill="none")+
          xlab(paste("Risk of",name_i,Unit_variable[[x[1]]]))+
          ylab(paste("Risk of",name_j,Unit_variable[[x[2]]]))+
          theme(axis.title=element_text(size=15),
                legend.text=element_text(size=10))
        ggsave(filename =paste0(Lg_path,"/Distrib_RiskF",name_i,"_",name_j,".png"),
               plot = GG_RF_j,width=8,height=6)
        dev.off()
})
F_automatic_ggplot_marg<-function(vect_time_choice){
  obj<-Marg_2d_simul_vs_obs(list_simul = L_simul,
                       list_obs = L_obs,
                       vect_times = vect_time_choice,
                       l_name_2V =L_nameV ,
                       cols_ggplot=cols_)
  ggsave(filename = paste0(TEND_path
                           ,"/distrib_values_",L_nameV[1],"t=",
                           vect_time_choice[1],"X",L_nameV[2],"t=",
                           vect_time_choice[2],CPLMT,".png"),
         width=8,height = 6,
         plot = obj)
}

# Evaluate the new expression
list_pairs<-list("1"=c(1,19),"2"=c(19,37),
                  "3"=c(1,37),
                 "4"=c(13,19))
sapply(list_pairs,F_automatic_ggplot_marg)

# Use of the machine learning -----------------------------------------------------
H_Params<-c("radial","logit",500)
NB_times<-100
ALPHA_PROP<-0.10
K<-4
typeS<-"simple"
resultat<-matrix(NA,ncol=3,nrow=3)
Lperfs_Y<-Running_perfs_ML(Base_simul = Simuls,Base_data = Reals, 
                           hyp_param = H_Params, NB_times =NB_times,alpha_prop=ALPHA_PROP,K = K, 
                           type_sampling = typeS,title_ROC =variable,
                           NB_shown_ROC = 20)

resultat[1,]<-Lperfs_Y$Accuracy$qu_accuracy
resultat

# End of the code ---------------------------------------------------------
#########
