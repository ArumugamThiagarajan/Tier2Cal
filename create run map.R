# 1 - sites with multiple TrtID per StockID
# run them by StockID, no problem

# Processing of site predictors (sand, till, cinput, etc.) for - sites with multiple StockID per treatmentID
# Predictors and C stocks should be aggregated by treatment


run.map <- function(site_data_ipcct2,stocks.df,dsm_data){
## Aggregation of predictors by treatment
 test.t <- site_data_ipcct2 %>% group_by(TrtID_Final) %>%
   summarise(stocks=length(unique(site))) %>%
   mutate(multi_stocks=ifelse(stocks>1,TRUE,FALSE))
 
 test.s <- site_data_ipcct2 %>% group_by(site) %>%
   summarise(trt=length(unique(TrtID_Final))) %>%
   mutate(multi_treatments=ifelse(trt>1,TRUE,FALSE))

 multi.trt.sites <- test.s %>% filter(multi_treatments==TRUE) %>% pull(site)
 multi.trt.treatments <- site_data_ipcct2 %>% filter(site %in% multi.trt.sites) %>% pull(unique(TrtID_Final))
 
## Aggregation of C stocks by best method to obtain longer time series
#trt_on <- "E1_F25NS_CT_CR1"

for (i in 1:nrow(test.t)){
  trt_on <- test.t[i,]$TrtID_Final
  mt <- test.t[i,]$multi_stocks
  
  if (mt==FALSE & trt_on %in% multi.trt.treatments){ # One site, multiple treatments (consecutive treatments, often change in crop rotation over time)
    pre_sub <- subset(stocks.df,stocks.df$TrtID_Final==trt_on)
    stock_on <- pre_sub[1,]$site
    sub <- subset(stocks.df,stocks.df$site==stock_on)
    
    out <- sub %>% 
      mutate(TrtID_Final=unique(TrtID_Final) %>% sort %>% paste(collapse = " + ")) %>%
      group_by(site,year) %>%
      summarise(TrtID_Final=first(TrtID_Final),
        Exp_ID=first(experiment),
                SOC_mean=mean(modelled_SOC,na.rm=T),
                SOC_sd=sd(modelled_SOC,na.rm=T),
                n = n(),
                .groups="keep") 
    
  } 
  if (mt==FALSE & !(trt_on %in% multi.trt.treatments)){ # One site, one treatment, simplest situation
    sub <- subset(stocks.df,stocks.df$TrtID_Final==trt_on)
    out <- sub %>% 
      group_by(site,year) %>%
      summarise(TrtID_Final=unique(TrtID_Final) %>% sort %>% paste(collapse = " + "),
                Exp_ID=first(experiment),
                SOC_mean=mean(modelled_SOC,na.rm=T),
                SOC_sd=sd(modelled_SOC,na.rm=T),
                n = n(),
                .groups="keep") ## Implement an equivalent for single-stock treatments
    
  } 
  if (mt==TRUE) { #Multiple sites with the same treatments
    sub <- subset(stocks.df,stocks.df$TrtID_Final==trt_on)
    out <- sub %>% 
      group_by(TrtID_Final,year) %>%
      summarise(site=unique(site) %>% sort %>% paste(collapse = " + "),
                Exp_ID=first(experiment),
                SOC_mean=mean(modelled_SOC,na.rm=T),
                SOC_sd=sd(modelled_SOC,na.rm=T),
                n = n(),
                .groups="keep")
  }
  if (!exists("stocks_agg_temp")){stocks_agg_temp <- out} else {stocks_agg_temp <- rbind(stocks_agg_temp,out)}
  rm(out)
} 

stocks_agg_temp <- stocks_agg_temp %>% distinct()
  
#Determine if time series will be run by treatment or by "UniqueStockID" (for consecutive treatments in the same experiment)

run.map <- stocks_agg_temp %>%
  ungroup() %>%
  select(site,TrtID_Final,Exp_ID) %>%
  rename(UniqueStockID=site) %>%
  distinct() %>%
  mutate(RunBy = ifelse(grepl("+",TrtID_Final,fixed = TRUE),
                        "stock",
                        "treatment"),
         till=NA,
         first.year=NA,
         initial.SOC=NA)

## Initial SOC 

for (i in 1:nrow(run.map)){
  run.by <- run.map[i,]$RunBy
  if (run.by=="treatment"){
    trt.on<-run.map[i,]$TrtID_Final
    sub.pred<-site_data_ipcct2 %>% filter(TrtID_Final==trt.on)
    sub.stock<-stocks.df %>% filter(TrtID_Final==trt.on)
    first.year<- min(sub.pred$year)
    initial.soc<- sub.stock[sub.stock$year==first.year,]$modelled_SOC
  }
  
  if (run.by=="stock"){
    stock.on<-run.map[i,]$UniqueStockID
    sub.pred<-site_data_ipcct2 %>% filter(site==stock.on)
    sub.stock<-stocks.df %>% filter(site==stock.on)
    first.year<- min(sub.pred$year)
    initial.soc<- sub.stock[sub.stock$year==first.year,]$modelled_SOC
  }
  run.map[i,]$till=first(sub.pred$till)
  run.map[i,]$first.year=first.year
  run.map[i,]$initial.SOC=ifelse(length(initial.soc)==0,NA,initial.soc)
  rm(first.year)
  rm(initial.soc)
  
}

# Add DSM data and impute missing initial C values (most) using the DSM data
run.map <- run.map %>% merge(dsm_data %>% select(Exp_ID,location_name,Can_SOC30),by="Exp_ID")
# Note: should label the initial soc clusters, they will become iterative parameters later on
run.map <- run.map %>% mutate(initial.SOC.impute= ifelse(!is.na(initial.SOC),initial.SOC,Can_SOC30))

return(run.map)
}





