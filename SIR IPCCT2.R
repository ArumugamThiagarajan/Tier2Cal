## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

#=================================================================================
# Sampling Importance Resampling algorithm with Latin hypercube sampling in R
#=================================================================================

run.SIR.IPCCT2 <- function(main.dir,
                           Experiment_name,
                           parameter.df,
                           stocks.df,
                           calibration.df,
                           climate_data,
                           RunMap,
                           sample_size,
                           resample_size) {
  Lkhood <- NULL
  Lkhood_list <- list()
  
  loglik=function(m,o){
    if(length(m)!=length(o)){
      print("Inequal number of modeled and observed values, cannot proceed")
      return()
    }
    
    res=log(m)-log(o)
    sigma=sqrt(mean(res^2))
    n=length(m)
    lk=-n*log(sigma)-(1/(2*sigma^2))*sum(res^2)
    return(lk)
    
  }
  
  RunIndex.list=unique(calibration.df$RunIndex)
  
  n_clust=detectCores()
  cl <- makeCluster(n_clust-4)
  registerDoParallel(cl)
  
  df_par <- foreach(RunIndex=RunIndex.list, .combine=rbind,  .inorder=TRUE) %dopar% {
    
    library(caret)
    library(dplyr)
    library(boot)
    library(doParallel)
    library(parallel)
    library(foreach)
    library(purrr)
    library(tidyverse)
    
    source(paste0(main.dir,"//scripts/Recent/ipcct2_SOCTier2Model.r"))
    source(paste0(main.dir,"//scripts//Recent//Spinnup IPCCT2.R")) # Model initiation for IPCCt2
    
    loglik=function(m,o){
      if(length(m)!=length(o)){
        print("Inequal number of modeled and observed values, cannot proceed")
        return()
      }
      
      res=log(m)-log(o)
      sigma=sqrt(mean(res^2))
      n=length(m)
      lk=-n*log(sigma)-(1/(2*sigma^2))*sum(res^2)
      return(lk)
      
    }
    
    for (i in 1:nrow(parameter.df)){
      
      parameters <- parameter.df[i,]
      
      RunIndex.sub=calibration.df[calibration.df$RunIndex==RunIndex,] %>%
        arrange(year) %>%
        distinct() %>%
        mutate(Diff = year -lag(year)) %>%
        replace_na(list(Diff=0))
      
      ## Address discontinuity in data (mostly due to crop rotation) ##!! Something to improve here !!##
      if (any(RunIndex.sub$Diff>1)) {
        problems=RunIndex.sub[RunIndex.sub$Diff>1,]
        RunIndex.sub <- RunIndex.sub %>% 
          slice(1:min(as.numeric(row.names(problems))-1))
      }
      
      spinup.df <- initiate.ipcc2(index.df=RunMap[RunMap$RunIndex==RunIndex,],
                                  params=parameters,
                                  site_data_ipcct2 = site_data_ipcct2,
                                  climate_data_ipcct2= climate_data_ipcct2)
      
      init_active = spinup.df[1,]$init_active
      init_slow = spinup.df[1,]$init_slow
      init_passive = spinup.df[1,]$init_passive
      
      climate_in <- climate_data[climate_data$POLYID==init.df [1,]$POLYID,]
      
      climate_normal <- climate_in %>% group_by(month) %>% summarise_all(mean) %>% select(month,tavg,mappet,irrig)
      
      if (min(RunIndex.sub$year)<1981){
        climate_int <- data.frame(POLYID=RunIndex.sub[1,]$POLYID,
                                  year=rep(min(RunIndex.sub$year):1980, each=12),
                                  month= rep(1:12,1981-min(RunIndex.sub$year))) %>%
          merge(climate_normal,by=("month"))
        climate_in <- rbind(climate_int,climate_in)
      }
      
      
      
      modelled <- IPCCTier2SOMmodel(SiteData= RunIndex.sub,
                                    wth = climate_in,
                                    init.active = init_active,
                                    init.slow = init_slow,
                                    init.passive = init_passive,
                                    params=parameters)
      
      if (init.df [1,]$RunBy=="treatment"){
        actuals <-  stocks.df %>%
          select(TrtID_Final, year,  actual = modelled_SOC)
        
        model_actual <- modelled %>%
          merge(actuals, by=c("TrtID_Final", "year")) %>%
          select(RunIndex, year, soc_total, actual) %>%
          filter(!is.na(actual))
        
        if (!is.na(init.df[1,]$initial.SOC)){model_actual <- model_actual %>% filter(year!=init.df [1,]$first.year)} #Remove the first year when it was measured (would provide perfect fit because model was initiated with this value)
        
      } 
      if (init.df [1,]$RunBy=="stock"){
        actuals <-  stocks.df %>%
          select(site, year,  actual = modelled_SOC)
        
        model_actual <- modelled %>%
          merge(actuals, by=c("site", "year")) %>%
          select(RunIndex, year, soc_total, actual) %>%
          filter(!is.na(actual))
        
        if (!is.na(init.df[1,]$initial.SOC)){model_actual <- model_actual %>% filter(year!=init.df [1,]$first.year)} #Remove the first year when it was measured (would provide perfect fit because model was initiated with this value)
        
      }
      
      loglike <- loglik(model_actual$soc_total, model_actual$actual)
      
      out.df=data.frame(parameters,
                        RunIndex=RunIndex,
                        loglike=loglike)
      
      if (i==1){final.df <- out.df} else {final.df=rbind(final.df,out.df)}
    }
    write.csv(final.df,paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//SIR//interim//",RunIndex,".csv"),row.names<-FALSE)
    return(final.df)
  }
  
  stopCluster(cl)
  Lkhood1 <- df_par %>%
    mutate(loglike = ifelse(loglike == -Inf, NA, loglike)) %>%
    group_by(SampleID) %>%
    summarise(loglike = mean(loglike, na.rm=T)) %>%
    mutate(weights = exp(loglike)/sum(exp(loglike)))
  
  sampIndx <- sample(1:nrow(Lkhood1), 
                     size = SIR_resample_size,
                     replace = FALSE,
                     prob = Lkhood1$weights)
  
  PostTheta <- as.data.frame(parameter.df[sampIndx,])
  return(PostTheta)
}


