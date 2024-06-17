## Authors: 
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023
#site_data <- read.csv(paste0(main.dir,"//data//input data//site data//LTE_Master_beta_sept_25.csv"))
#dsm_data<- read.csv(paste0(main.dir,"//data//input data//site data//lte_soc_30_cm_dsm_expid_may_2.csv"))
#merge.key="location_name"
get.stocks <- function(site_data) {
  stocks.df <- site_data %>%
    mutate(
      site = UniqueStockID,
      year = year_name,
      depthID = paste0(soil_depth_min_cm,"-",soil_depth_max_cm)) %>%
    select(site,
           year,
           Exp_ID,
           location_name,
           TrtID_Final,
           soc_tha,
           soil_total_carbon_px,	
           soil_total_carbon2_px,
           bd_gcm3,
           soil_depth_min_cm,
           soil_depth_max_cm,
           depthID) %>%
    distinct() %>%
    arrange(site,year) %>%
    mutate(soc_tha=as.numeric(soc_tha))
  
  # Get C content of soil layer from soc_tha or from bulk density x carbon percent
  
  stocks.df <- stocks.df %>%
    mutate(soc_tha=ifelse(soc_tha==0,NA,soc_tha))%>%
    mutate(c_percent_mean=ifelse(is.na(soil_total_carbon_px)&is.na(soil_total_carbon2_px),
                                 NA,
                                 (soil_total_carbon_px+soil_total_carbon2_px)/2)) %>%
    mutate(c_layer=case_when(
                             !is.na(soc_tha)~soc_tha,
                             !is.na(c_percent_mean) & !is.na(bd_gcm3) ~ c_percent_mean/100*bd_gcm3*(soil_depth_max_cm-soil_depth_min_cm)*100))
  
  missing.SOC.list <- stocks.df %>% 
    group_by(TrtID_Final) %>% 
    summarise(mean_soc=mean(c_layer,na.rm=T)) %>% 
    filter(is.nan(mean_soc)) %>% 
    select(TrtID_Final)
  if(nrow(missing.SOC.list)>0){print(paste0("no C stock measurements available for treatment: ",missing.SOC.list$TrtID_Final))}
  
  # remove NAs
  stocks.df <- stocks.df %>%
    filter(!is.na(c_layer))

  return (stocks.df)
}


obtain.dsm.correction <- function(stocks.df,dsm_data,merge.key) {

# Join the DSM data to the stocks
dsm=merge(stocks.df,dsm_data,by=merge.key,left=T) %>%
  select(location_name,Can_SOC30) %>% 
  group_by(location_name) %>%
  summarise(Can_SOC30=mean(Can_SOC30))

## Get correction factors for top 30 cm for each Exp_ID and create QA flag
## These factor will not necessarily be used
correct.fac <- stocks.df %>%
  merge(dsm,by=merge.key,left=T) %>%
  mutate(dif=c_layer-Can_SOC30,min=soil_depth_min_cm,max=soil_depth_max_cm) %>%
  group_by(location_name,depthID) %>%
  arrange(year) %>%
  slice(1) %>%
  group_by(depthID,location_name,min,max) %>%
  summarise(SOC_correction_mean=mean(dif,na.rm=T)*-1,
            SOC_correction_median=median(dif,na.rm=T)*-1,
            .groups="keep") %>%
  filter(!is.na(SOC_correction_mean)) %>%
  mutate(QA_flag1=case_when(min==0 & max<30 & SOC_correction_mean>0 ~ 1,
                            min==0 & max<30 & SOC_correction_mean<0 ~ 0,
                            min==0 & max>30 & SOC_correction_mean<0 ~ 1,
                            min==0 & max>30 & SOC_correction_mean>0 ~ 0)) 

return(correct.fac)
}


#test 12-12
# Check the counts of depth measurements
#stocks.df$depthID=paste0(stocks.df$soil_depth_min_cm,"-",stocks.df$soil_depth_max_cm)
#table(stocks.df$depthID)

# Check the counts of depth measurements by experiment
#table(stocks.df$depthID,stocks.df$Exp_ID)

# For debugging...
#df=stocks.df
#site=unique(df$site)[63]#
#site=unique(df$site)[211] #
#site="E7_Field 4T9"
#year=1994
#year=2003
#depth_cm=30



# Function to estimate 30 cm SOC at all sites and all dates
top.SOC.estimate <- function(df,correct.fac,depth_cm){
  df <- df[!is.na(df$c_layer),]
  df <- df[df$c_layer!=0,]
  
  for (site in unique(df$site)){
    # select one site
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      correction<- 0 #(re)-initialise correction factor
      # Select one year, get all available measurements for that year
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ") # index name of the depth increment
      TrtID_Final=unique(sub.df2$TrtID_Final)
      
      #Create empty data frame with rows for each 0.5 cm increment
      res.df=data.frame(depth.min=seq(0,depth_cm-0.5,0.5),
                        depth.max=seq(0.5,depth_cm,0.5),
                        measured_c=NA)
      
      
      
      # Adapt the 
      if (sub.df2[1,]$soil_depth_min_cm==0 & sub.df2[1,]$soil_depth_max_cm>depth_cm | # Catch cases with measurements exceeding the target depth increment (e.g., 0-60 when we want 0-30)
          max(sub.df2$soil_depth_max_cm)<depth_cm) { # Catch cases with measurements short of the target depth increment (e.g., 0-20 when we want 0-30)
        
        cor.sub=correct.fac %>% filter(location_name==sub.df2[1,]$location_name, min==min(sub.df2$soil_depth_min_cm), max==max(sub.df2$soil_depth_max_cm))
        if (is.na(cor.sub[1,]$QA_flag1)) { #Check the quality of the correction factor
          correction<-0 #Don't apply any correction if it we can't trust the estimated factor
        } else if (cor.sub[1,]$QA_flag1==1) {correction<-cor.sub$SOC_correction_mean}
      }
      if (max(sub.df2$soil_depth_max_cm)>depth_cm){ #process for long
        modelled_SOC=sum(sub.df2[1,]$c_layer, # First row will have the soc total
                         correction, # apply depth correction from dsm estimate
                         na.rm=T)
      } else {
      # for all other measurements (correct or short), fill the table
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$c_layer/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)}
      
    res.df <- replace(res.df, res.df==0, NA)
    
    res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T) # Take the row mean for each increment
    
    modelled_SOC=sum(res.df$measured_c, # Get the sum of all measurements
                     correction, # apply depth correction from dsm estimate
                     na.rm=T)  
  
      }
    out=data.frame(site,
                   experiment=sub("\\_.*", "", site),
                   TrtID_Final,
                   year,
                   input_type=type_string,
                   modelled_SOC)
    
    if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
  }}
  return(final)
}

#site<-"E4_102"
#year<-1954
## Older function for comparison
top.SOC.raw <- function(df){
  df <- df[!is.na(df$c_layer),]
  df <- df[df$c_layer!=0,]
  
  for (site in unique(df$site)){
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ")
      TrtID_Final=unique(sub.df2$TrtID_Final)
      if (length(TrtID_Final)>1){print(paste0("error: more than one treatment ID for stockID --> ",site))}
      d.max=max(sub.df2$soil_depth_max_cm)
      res.df=data.frame(depth.min=seq(0,d.max-0.5,0.5),
                        depth.max=seq(0.5,d.max,0.5),
                        measured_c=NA)
      
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$c_layer/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)
      }
        res.df <- replace(res.df, res.df==0, NA)
        
        res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T)
        
        modelled_SOC=sum(res.df$measured_c)
        
      out=data.frame(site,
                     experiment=sub("\\_.*", "", site),
                     TrtID_Final,
                     year,
                     input_type=type_string,
                     modelled_SOC)
      
      if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
      
    }
  }
  return(final)
}

top.SOC.naive <- function(df,depth_cm){
  df <- df[!is.na(df$c_layer),]
  df <- df[df$c_layer!=0,]
  
  for (site in unique(df$site)){
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ")
      TrtID_Final=unique(sub.df2$TrtID_Final)
      if (length(TrtID_Final)>1){print(paste0("error: more than one treatment ID for stockID --> ",site))}
      d.max=max(sub.df2$soil_depth_max_cm)
      res.df=data.frame(depth.min=seq(0,depth_cm-0.5,0.5),
                        depth.max=seq(0.5,depth_cm,0.5),
                        measured_c=NA)
      
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$c_layer/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)
        
        res.df <- replace(res.df, res.df==0, NA)
        
        res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T)
        
        # Filling missing values in top / bottom section(s)
        res.df <- res.df %>%
        fill(measured_c, .direction = "down")%>%
        fill(measured_c, .direction = "up")
        
        modelled_SOC=sum(res.df$measured_c)
      }
      out=data.frame(site,
                     experiment=sub("\\_.*", "", site),
                     TrtID_Final,
                     year,
                     input_type=type_string,
                     modelled_SOC)
      
      if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
      
    }
  }
  return(final)
}
