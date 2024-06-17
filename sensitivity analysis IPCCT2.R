## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

library(caret)
library(dplyr)
library(sensitivity)
library(boot)
library(doParallel)
library(parallel)
library(foreach)
library(purrr)
library(tidyverse)
library(ggplot2)

run.GSA.IPCCT2 <- function(main.dir,
                           parameter_bounds,
                           stocks_data,
                           calibration_data,
                           site_data_ipcct2,
                           RunMap,
                           climate_data_ipcct2,
                           method,
                           sample_size,
                           likelihood="stock") {
  # read prior distribution from a csv file
  # (Required columns: Parameter, value, lower, upper)
  paramBounds <- parameter_bounds %>%
    arrange(Parameter)
  
  # names of parameters that are allowed to vary
  varSI       <- paramBounds$Parameter
  nParams     <- length(varSI)
  
  # sample size (10 used for illustration purposes)
  # (1024 used in Gurung et al., 2020)
  N <- sample_size
  
  # Sobols method required 2 random matrix
  m1 = matrix(runif(nParams*N), nrow=N);
  m2 = matrix(runif(nParams*N), nrow=N);
  M1 <- matrix(0, nrow = N, ncol = nParams)
  M2 <- matrix(0, nrow = N, ncol = nParams)
  
  # transform standard uniform to prior distribution 
  for(i in 1:nParams){
    pos <- which(paramBounds$Parameter == varSI[i])
    lower <- paramBounds[pos, "lower"]
    upper <- paramBounds[pos, "upper"]
    M1[, i] <- qunif(m1[, i], min = lower, max = upper)
    M2[, i] <- qunif(m2[, i], min = lower, max = upper)
  }
  X1 = data.frame(M1)
  X2 = data.frame(M2)
  names(X1) <- varSI
  names(X2) <- varSI
  
  if(method == "fast99") {
    fast99_qargs <- paramBounds %>%
      group_by(Parameter) %>%
      select(Parameter, min = lower, max = upper) %>%
      group_split %>%
      map(function(x) {y <- x %>%
        select(-Parameter) %>%
        as.list()
      })
  }
  
  if (method=="soboljansen"){
    si_obj2 <- sensitivity::soboljansen(model = NULL, X1 = X1, X2 = X2, nboot = 100)
  }
  if (method=="fast99"){
    si_obj2 <- sensitivity::fast99(model = NULL,
                                                   factors = paramBounds$Parameter,
                                                   n = N, M = 4,
                                                   q="qunif",
                                                   q.arg = fast99_qargs)
  }
  
  X <- si_obj2$X
  X <- cbind("SampleID" = 1:nrow(X), X)
  # NaN values can be in X if sample size was too small.
  if(any(is.nan(unlist(X)))) {warning("gsa: NaN values detected in sensitivity analysis. Try increasing sample_size")}
  
  
  params_list_sorted_names <- c("SampleID",varSI)
  params_list <- X %>%
    rowwise %>%
    group_split %>%
    map(function(x) {
      y <- x %>%
        pivot_longer(everything()) %>%
        deframe()
      z <- split(unname(y),names(y))
      return(z)
    }) %>%
    map(~ .[params_list_sorted_names])
  
  RunIndex.list=unique(calibration_data$RunIndex)
  
  allConnections<-showConnections(all = FALSE)
  
  if (length(allConnections) > 0) {
    # Close all open connections
    closeAllConnections()
    cat("Closed all open connections before starting....\n")
  }
  
  n_clust=detectCores()
  cl <- makeCluster(n_clust-4)
  registerDoParallel(cl)
  
  #Debug RunIndex=RunIndex.list[43]
  df_par <- foreach(RunIndex=RunIndex.list, .combine=rbind,  .inorder=TRUE) %dopar% {
    
    library(caret)
    library(dplyr)
    library(sensitivity)
    library(boot)
    library(doParallel)
    library(parallel)
    library(foreach)
    library(purrr)
    library(tidyverse)
    library(imputeTS)
    
    source(paste0(main.dir,"//scripts//Recent//ipcct2_SOCTier2Model.r"))
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
    
  for (i in 1:nrow(X)){
    RunIndex.sub=calibration_data[calibration_data$RunIndex==RunIndex,] %>%
      arrange(year) %>%
      distinct() %>%
      mutate(Diff = year -lag(year)) %>%
      replace_na(list(Diff=0))
    
    parameters = X[i,]

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
                                  merge(climate_normal,by="month")
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
    } else if (likelihood=="delta"){
    loglike <- loglik.delta(m=model_actual$delta_soc_model,
                            o=model_actual$delta_soc_actual)  
    } else if (likelihood=="cumul") {
    loglike <- loglik.delta(m=model_actual[nrow(model_actual),]$cumul_soc_model,
                            o=model_actual[nrow(model_actual),]$cumul_soc_actual)
    }
    
    out.df=data.frame(parameters,
                      RunIndex=RunIndex,
                      loglike=loglike)
    
    
    if (i==1){final.df <- out.df} else {final.df=rbind(final.df,out.df)}
    }
    return(final.df)
  }
  
  stopCluster(cl)
  
  Lkhood1 <- df_par %>%
    mutate(loglik = ifelse(loglike == -Inf, NA, loglike)) %>%
    group_by(id=SampleID) %>%
    summarise(loglik = mean(loglik, na.rm=T))
  
  si_obj2_llkhd <- sensitivity::tell(x = si_obj2, y = Lkhood1$loglik)
  
  
  if(method == "soboljansen"){
    # Calculate First-order and Total global sensitivity indices
    singleSI <- si_obj2_llkhd$S %>%
      select(SI = original, lci = `min. c.i.`, uci = `max. c.i.`) %>%
      rownames_to_column("params") %>%
      mutate(type="single")
    
    totalSI <- si_obj2_llkhd$T %>%
      select(SI = original, lci = `min. c.i.`, uci = `max. c.i.`) %>%
      rownames_to_column("params") %>%
      mutate(type="Total")
    
    combined_si <- rbind(singleSI,totalSI)
    
    return_si <- data.frame(combined_si)
  }
  
  if(method == "fast99"){
    return_si <- tibble(main = si_obj2_llkhd$D1 / si_obj2_llkhd$V,
                        interactions = 1 - si_obj2_llkhd$Dt / si_obj2_llkhd$V) %>%
      mutate(params =  varSI, .before=main)
  }
  
  return(return_si)
}



