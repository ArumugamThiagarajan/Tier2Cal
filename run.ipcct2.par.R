
run.IPCCT2.par <- function(main.dir,
                           results.dir,
                           Experiment_name,
                           out.subfolder,
                           parameter.df,
                           stocks_data,
                           calibration.df,
                           site_data_ipcct2,
                           climate_data_ipcct2,
                           RunMap,
                           out_type="fit",
                           likelihood="stock") {
  
  RunIndex.list=unique(calibration.df$RunIndex)
  
  allConnections<-showConnections(all = FALSE)
  
  if (length(allConnections) > 0) {
    # Close all open connections
    closeAllConnections()
    cat("Closed all open connections before starting....\n")
  }
  
  n_clust=detectCores()
  cl <- makeCluster(n_clust-4)
  registerDoParallel(cl)
  
  #RunIndex=RunIndex.list[1]
foreach(RunIndex=RunIndex.list) %dopar% {
    
  
  library(caret)
  library(dplyr)
  library(boot)
  library(doParallel)
  library(parallel)
  library(foreach)
  library(purrr)
  library(tidyverse)
  library(imputeTS)
  
  source(paste0(main.dir,"//scripts/Recent/ipcct2_SOCTier2Model.r"))
  source(paste0(main.dir,"//scripts//Recent//Spinnup IPCCT2.R")) # Model initiation for IPCCt2
    
  loglik.stock=function(m,o){
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
  
  
  loglik.delta=function(m,o){
    if(length(m)!=length(o)){
      print("Inequal number of modeled and observed values, cannot proceed")
      return()
    }
    
    n = length(m)
    residual = m - o
    sigma = sqrt(sum(residual^2) / (n+1 - 1))  # Estimate sigma
    lk = -n * log(2 * pi * sigma^2) / 2 - sum((m - o)^2) / (2 * sigma^2)
    
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
      #if (any(RunIndex.sub$Diff>1)) {
        #problems=RunIndex.sub[RunIndex.sub$Diff>1,]
        #RunIndex.sub <- RunIndex.sub %>% 
          #slice(1:min(as.numeric(row.names(problems))-1))
      #}
      
      spinup.df <- initiate.ipcc2(index.df=RunMap[RunMap$RunIndex==RunIndex,],
                                  params=parameters,
                                  site_data_ipcct2 = site_data_ipcct2,
                                  climate_data_ipcct2= climate_data_ipcct2)
      
      init_active = spinup.df[1,]$init_active
      init_slow = spinup.df[1,]$init_slow
      init_passive = spinup.df[1,]$init_passive
      
      init.df <- RunMap[RunMap$RunIndex==RunIndex,]
      
      climate_in <- climate_data_ipcct2[climate_data_ipcct2$POLYID==init.df [1,]$POLYID,]
      
      climate_normal <- climate_in %>% group_by(month) %>% summarise_all(mean) %>% select(month,tavg,mappet,irrig)
      
      if (min(RunIndex.sub$year)<1981){
        climate_int <- data.frame(POLYID=RunIndex.sub[1,]$POLYID,
                                  year=rep(min(RunIndex.sub$year):1980, each=12),
                                  month= rep(1:12,1981-min(RunIndex.sub$year))) %>%
          merge(climate_normal,by=("month"))
        climate_in <- rbind(climate_int,climate_in)
      }
      
      modelled <- IPCCTier2SOMmodel(SiteData= RunIndex.sub %>% rename(site=UniqueStockID,cinput=cinput_t),
                                    wth = climate_in,
                                    init.active = init_active,
                                    init.slow = init_slow,
                                    init.passive = init_passive,
                                    params=parameters)
      
      if (init.df [1,]$RunBy=="treatment"){
        actuals <-  stocks_data %>%
          select(TrtID_Final, year,  actual = modelled_SOC)%>%
          filter(TrtID_Final==init.df[1,]$TrtID_Final)%>%
          add_row(TrtID_Final=init.df[1,]$TrtID_Final,
                  year=init.df[1,]$first.year,
                  actual=init.df[1,]$init.soc.mu,
                  .before = 1)
        
        model_actual <- modelled %>%
          merge(actuals, by=c("TrtID_Final", "year")) %>%
          select(RunIndex, year, soc_total, actual) %>%
          filter(!is.na(actual))
      } else if (init.df [1,]$RunBy=="stock"){
        actuals <-  stocks_data %>%
          select(site, year,  actual = modelled_SOC)%>%
          filter(site==init.df[1,]$UniqueStockID)%>%
          add_row(site=init.df[1,]$UniqueStockID,
                  year=init.df[1,]$first.year,
                  actual=init.df[1,]$init.soc.mu,
                  .before = 1)
        
        model_actual <- modelled %>%
          merge(actuals, by=c("site", "year")) %>%
          select(RunIndex, year, soc_total, actual) %>%
          filter(!is.na(actual))
      }
      
      model_actual <- model_actual %>%
        distinct() %>%
        mutate(lag_soc_actual=lag(actual),
               lag_soc_model=lag(soc_total)) %>%
        mutate(delta_soc_actual=actual-lag_soc_actual,
               delta_soc_model=soc_total-lag_soc_model) %>%
        mutate(delta_soc_actual=na_replace(delta_soc_actual,0),
               delta_soc_model=na_replace(delta_soc_model,0)) %>%      
        mutate(cumul_soc_actual=cumsum(delta_soc_actual),
               cumul_soc_model=cumsum(delta_soc_model))
      
      
      model_actual <- model_actual %>% filter(year!=init.df [1,]$first.year) #Remove the first year (would provide perfect fit because model was initiated with this value)
      
      if (likelihood=="stock"){
        loglike <- loglik.stock(m=model_actual$soc_total,
                                o=model_actual$actual)
        RMSE <- Metrics::rmse(actual=model_actual$actual,
                              predicted= model_actual$soc_total)
        bias <- Metrics::bias(actual=model_actual$actual,
                              predicted= model_actual$soc_total)
        
      } else if (likelihood=="delta"){
        loglike <- loglik.delta(m=model_actual$delta_soc_model,
                                o=model_actual$delta_soc_actual)  
        RMSE <- Metrics::rmse(actual=model_actual$delta_soc_actual,
                              predicted= model_actual$delta_soc_model)
        bias <- Metrics::bias(actual=model_actual$actual,
                              predicted= model_actual$soc_total)
      } else if (likelihood=="cumul") {
        loglike <- loglik.delta(m=model_actual[nrow(model_actual),]$cumul_soc_model,
                                o=model_actual[nrow(model_actual),]$cumul_soc_actual)
        RMSE <- Metrics::rmse(actual=model_actual[nrow(model_actual),]$cumul_soc_actual,
                              predicted= model_actual[nrow(model_actual),]$cumul_soc_model)
        bias <- Metrics::bias(actual=model_actual$actual,
                              predicted= model_actual$soc_total)
      }
      
      out.df=data.frame(parameters,
                        RunIndex=RunIndex,
                        loglike=loglike,
                        RMSE=RMSE,
                        bias=bias)
      
      if (i==1){final.df <- out.df} else {final.df=rbind(final.df,out.df)}
      if (i==1){final.model_actual <- model_actual} else {final.model_actual=rbind(final.model_actual,model_actual)}
    }
    if (out_type=="fit"){
    write.csv(final.df,paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//",out.subfolder,RunIndex,".csv"),row.names<-FALSE)
    }
    if (out_type=="values"){
      write.csv(final.model_actual,paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//",out.subfolder,RunIndex,".csv"),row.names<-FALSE)
    }
      }
stopCluster(cl)
}
