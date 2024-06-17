data.prep.ipcct2 <- function (site_data_raw,
                              cinput.df,
                              poly.df) {
  
library(imputeTS)
  
site_data_ipcct2 <- site_data_raw %>%
  merge(poly.df, by=c("Exp_ID")) %>%
  mutate(tillage = ifelse(is.na(tillage), "CT", tillage)) %>%
  mutate(tillage = ifelse(tillage == "CT", "FT", tillage)) %>%
  mutate(
    POLYID = POLYID,
    year = year_name,
    sand = sand_px/100,
    clay = clay_px/100,
    till = tillage,
    irrig=0) %>%
  rowwise() %>% 
  select(UniqueStockID,Exp_ID,POLYID,TrtID_Final,year,sand,clay,till,irrig,AAFC_code) %>%
  group_by(UniqueStockID,POLYID,Exp_ID,TrtID_Final,year,till,AAFC_code) %>% #Get mean values when multiple depth measurements are available
  summarise(sand=mean(sand,na.rm=T),
            clay=mean(clay,na.rm=T),
          irrig=mean(irrig,na.rm=T),
          .groups="keep")%>% 
  mutate(AAFC_code_adapted=AAFC_code) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+BROME","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+CLOV","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+CWG+OHAYFD","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="BARLEY+ALFALFA+BROME","BARLEY",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="BARLEY+CLOV","BARLEY",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="FALLOW","SUMMRF",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="SUMMRFG","SUMMRF",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="CSFG","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="GML","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="GRS","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="OAT+OHAYFD","OATS",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="OATS+PEA","OATS",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="PEA","PEAS",AAFC_code_adapted))

# Remove the treatments containing perennial crops
perenials=c("CLOVER","CLOV","CWG+HAY","ALFALFA","OHAYFD","GRS")

site_data_ipcct2 <- site_data_ipcct2 %>% group_by(TrtID_Final) %>%
  mutate(contains_peren = ifelse(any(AAFC_code_adapted %in% perenials),TRUE,FALSE)) %>%
  ungroup %>%
  filter(contains_peren ==FALSE) %>%
  select(-contains_peren)

site_data_ipcct2 <- site_data_ipcct2 %>% left_join(cinput.df,by=c("TrtID_Final","UniqueStockID","year"))

interpolated.df<- foreach(i=unique(site_data_ipcct2$UniqueStockID),.combine='rbind')%do%{
  sub.df=site_data_ipcct2[site_data_ipcct2$UniqueStockID==i,]
  int.df=data.frame(year=seq(min(sub.df$year),max(sub.df$year),1))
  int.df<- left_join(int.df,sub.df,by="year",multiple="first") 
  int.df2 <- int.df %>%
    mutate(UniqueStockID = first(UniqueStockID),
           POLYID = first(POLYID),
           Exp_ID = first(Exp_ID)) %>%
    fill(till,.direction="down") %>%
    fill(sand,.direction="down") %>%
    fill(clay,.direction="down") %>%
    fill(TrtID_Final,.direction="down") %>%
    mutate(irrig = na_interpolation(irrig)) 
  return(int.df2)
}

## Interpolation of missing data

interpolate.cinput <- function(x) {
  if (length(na.omit(x))>=10) {out=na_seasplit(x,find_frequency=TRUE)} 
  else if (length(na.omit(x))>=3) {out=na_kalman(x)} 
  else {out=na_replace(x,fill=mean(x,na.rm=T))}
  return(out)
}


interpolated.df <- interpolated.df %>%
  group_by(UniqueStockID) %>%
  arrange(year) %>%
  mutate(nfrac = interpolate.cinput(nfrac),
         ligfrac =interpolate.cinput(ligfrac),
        AGR_t= interpolate.cinput(AGR_t),
         BGR_t = interpolate.cinput(BGR_t),
         cinput_t = interpolate.cinput(cinput_t)) %>%
  ungroup()


return (interpolated.df)
}

################# Climate data ###################

data.prep.climate <- function (dir) {
  
  climate.dir <- dir
  # Climate data
  ## open and combine files
  files <- list.files(climate.dir, full.names = TRUE, pattern = "csv$")
  
  climate_data <- files %>%
    lapply(read.csv) %>% 
    bind_rows %>%
    rename(Tavg = Tmean,
           PREC = Precip,
           JulianDay = Julian)
  
  ## Add potential evapotransipration
  holos_calculate_pet <- function(meanDailyTemperature,
                                  solarRadiation,
                                  relativeHumidity){
    # This is a vectorized implementation of the Holos reference PET calculator.
    # https://github.com/holos-aafc/Holos/blob/main/H.Core/Calculators/Climate/EvapotranspirationCalculator.cs
    
    term1 = 0.013
    term2 = meanDailyTemperature / (meanDailyTemperature + 15)
    term3 = (23.8856 * solarRadiation) + 50
    term4 = 1 + ((50 - relativeHumidity) / 70)
    
    result <- ifelse(relativeHumidity >= 50,
                     term1 * term2 * term3,
                     term1 * term2 * term3 * term4)
    result <- ifelse(result < 0, 0, result)
    result <- ifelse(meanDailyTemperature <= 0, 0, result)
    
    return(result)
  }
  
  climate_data <- climate_data %>%
    mutate(PET = holos_calculate_pet(meanDailyTemperature=Tavg,
                                     solarRadiation=Rad,
                                     relativeHumidity=RH)) 
  
  # Format climate data and add irrigation to this dataset (required for IPCCT2 wth file)
  climate_data_ipcct2 <- climate_data %>%
    group_by(POLYID, year = Year, month = Month) %>%
    summarise(tavg = mean(Tavg),
              mappet = ifelse(is.infinite(sum(PREC) / sum(PET)),0,sum(PREC) / sum(PET))
    ) %>%
    mutate(irrig=0) # at this point we don't have any irrigation data, so set it to zero
  
  return(climate_data_ipcct2)
}
