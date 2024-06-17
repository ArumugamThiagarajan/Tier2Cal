
### Debug NP 03-24-23
#SiteData = site_data_ipcct2
#wth = climate_data_ipcct2
#init.active = init_active
#init.slow = init_slow
#init.passive = init_passive


# DEFAULT -----------------------------------------------------------------
IPCCTier2SOMmodel <- function(SiteData, wth, init.active,
                              init.slow,init.passive,
                              params) {
  
# Parameter overrides ---------------------------------------
if (length(params)==0){
 tillfac_FT = 3.036
 tillfac_RT = 2.075
 wfac_irri  = 0.775
 k10        = 18.5
 k20        = 4.9
 k30        = 7.4
 k40        = 0.209
 k50        = 0.00689
 f1         = 0.378
 f2         = 0.368
 f3         = 0.455
 f5         = 0.0855
 f6         = 0.0504
 f7         = 0.42
 f8         = 0.45
 tmax       = 45
 topt       = 33.69
 plig       = 3 } else {
   f1 = params[1,]$f1
   f2 = params[1,]$f2
   f3 = params[1,]$f3
   f5 = params[1,]$f5
   f6 = params[1,]$f6
   f7 = params[1,]$f7
   f8 = params[1,]$f8
   k10 = params[1,]$k10
   k20 = params[1,]$k20
   k30 = params[1,]$k30
   k40 = params[1,]$k40
   k50 = params[1,]$k50
   plig = params[1,]$plig
   tillfac_FT = params[1,]$tillfac_FT
   tillfac_RT = params[1,]$tillfac_RT
   tmax = params[1,]$tmax
   topt = params[1,]$topt
   wfac_irri = params[1,]$wfac_irri
 }
#browser()
  #  AUTHOR: 
  #    Ram Gurung and Stephen Ogle
  #    Natural Resource Ecology Laboratory
  #    Colorado State University
  #    Fort Collins, CO 80523
  #
  #    January 5, 2018
  #
  #  PURPOSE:
  #    function to estimate annual change in mineral SOC stocks
  #    based on a three pool steady state solution of the  
  #    Century model. (Reference: Paustian et al. 1997, Ogle et al. 2012)
  #
  #  INPUTS:
  #	   SiteDataFile	site level data.frame object with the following columns:
  #						site	- a unique identifier for the site (numerical or character)
  #						year	- simulation year
  #						sand	- sand fraction of the soil (0-1)
  #						cinput	- carbon input (g C/m^2)
  #						ligfrac	- lignin fraction of the cinput (0-1)
  #						nfrac	- nitrogen fraction of the cinput (0-1)
  #						till	- tillage practice one of the following
  #									FT: full-till 
  #									RT: reduced-till
  #									NT: no-till
  #	   WthFile weather file for the site with the following columns:
  #						year	- simulation year
  #						month	- simulation month (1-12)
  #						tave	- average monthly temperature (degree C))
  #						mappet	- average monthly precipitation divided by
  #								  average monthly PET
  #						irrig	- boolean (0/1) for irrigation for the month
  #									1 if irrigated
  #									0 if non-irrigated
  #	   ParamFile	a numerical vector of length 33 containing model parameters
  #					(see Line-57:86 for parameter descriptions)
  #	
  #  
  #  OUTPUTs are a data.frame object with the following columns
  #						site		- a unique identifier for the site as in the SiteData
  #						year		- simulation year
  #						dfac		- climate decomposition factor of the year for the site
  #						kA			- decay rate for active pool for the year
  #						SS_active	- steady state active pool (g C/m^2)
  #						active		- active pool (g C/m^2) for the year
  #						kS			- decay rate for slow pool for the year
  #						SS_slow		- steady state slow pool (g C/m^2)
  #						slow		- slow pool (g C/m^2) for the year
  #						kP			- decay rate for passive pool for the year
  #						SS_passive	- steady state passive pool (g C/m^2)
  #						passive		- passive pool (g C/m^2) for the year
  #						TotSOC		- total soil organic matter pool (g C/m^2)
  #									  (sum of active, slow, and passive pool for the year)
  #           deltaC - annual change in somsc (g C/m^2/year)
  #
  #=========================================================================================
  # MODEL PARAMETERS
  # tillfac_FT   = params$tillfac_ft               # tillage disturbance modifier tilled soil (Full Till)
  # tillfac_RT   = params$tillfac_rt               # tillage disturbance modifier tilled soil (Reduce Till)
  # wfac_irri    = params$wfacpar1               # wfac for irrigated field during the irrigation period
  # k10          = params$k10               # decay rate under optimum conditions for metabolic litter pool
  # k20          = params$k20             # decay rate under optimum condition for structural litter pool
  # k30          = params$k30               # decay rate under optimum condition for active
  # k40          = params$k40               # decay rate under optimum condition for slow
  # k50          = params$k50               # decay rate under optimum condition for passive
  # f1           = params$f1               # stabilization efficiencies for metabolic decay products entering the active pool 
  # f2           = params$f2              # stabilization efficiencies for structural decay products entering the active pool
  # f3           = params$f3              # stabilization efficiencies for structural decay products entering the slow pool
  # f5           = params$f5              # stabilization efficiencies for active pool decay products entering the passive pool
  # f6           = params$f6              # stabilization efficiencies for slow pool decay products entering the passive pool
  # f7           = params$f7              # stabilization efficiencies for slow pool decay products entering the active pool
  # f8           = params$f8              # stabilization efficiencies for passive pool decay products entering the active pool
  # tmax         = params$tmax              # maximum temperature on decomposition
  # topt         = params$topt              # optimum temperature on decomposition
  # plig         = params$plig              # empirical parameter to modify k20
  #=========================================================================================
  #  END OF MODEL PARAMETERS
  #=========================================================================================
  # Check for Consecutive Year
  years  <- sort(SiteData$year)
  ydiff      <- diff(years)
  yflag      <- sum(ydiff != 1)
  if(yflag != 0){
    #stop(paste("consecutive year needed for site ", SiteName, " and TreatmentID ", trt, sep = "")) #original
    stop(paste("consecutive year needed for site ", site_n, sep = "")) #modifed NP 03-24-23
  }
  # check for columns headings
  cflag1 <- sum(!is.element(c("site", "year", "sand", "cinput", "ligfrac", "nfrac", "irrig", "till"), names(SiteData)))
  if(cflag1 > 1){
    stop("check your SiteData file format.")
  }
  #=========================================================================================
  # Check completeness in weather data
  # check for years
  wflag <- sum(is.element(years, wth$year))
  if(wflag != length(years)){
    stop("wth file must have all the data associated with years in SiteData")
  }
  # check for columns headings
  cflag2 <- sum(!is.element(c("year", "month", "tavg", "mappet"), names(wth)))
  if(cflag2 > 1){
    stop("check your wth file format.")
  }
  #=========================================================================================
  # Check for parameter length
  # Note: This has been taken out since I changed the format of the input parameters to be wide -Francis
  # 
  #if(length(params$value) != 18){
    #stop("there must be 18 parameters.")
  #}
  #=========================================================================================
  # Other Check can be added here
  #=========================================================================================
  SOMstocks <- NULL
  # Initial SOM pool stocks
  Isom1   <- init.active
  Isom2   <- init.slow
  Isom3   <- init.passive
  # 
  for(year in years){
    #
    # cinput, management, and weather data for the year
    site       <- SiteData$site[SiteData$year == year]
    TrtID_Final       <- SiteData$TrtID_Final[SiteData$year == year]
    RunIndex <- SiteData$RunIndex[SiteData$year == year]
    cinput     <- SiteData$cinput[SiteData$year == year]
    TILL       <- SiteData$till[SiteData$year == year]
    sand       <- SiteData$sand[SiteData$year == year]
    L          <- SiteData$ligfrac[SiteData$year == year]    # lignin fraction of cinput
    N          <- SiteData$nfrac[SiteData$year == year]      # nitrogen fraction of cinput
    mtemp      <- as.numeric(wth$tavg[wth$year == year])
    mappet     <- as.numeric(wth$mappet[wth$year == year])
    IRRIG      <- as.numeric(wth$irrig[wth$year == year])
    #
    # estimate steady state soc pools and decay rates
    #=========================================================================================
    # Calculating AvgDfac (average decomposition factor )
    #=========================================================================================
    if(length(mappet) != 12){
      stop("Length of mappet must be 12. It is not 12, rather it is: ", length(mappet))
    }
    if(length(mtemp) != 12){
      stop("Length of monthly average temperature must be 12. It is not 12, rather it is: ", length(mtemp))
    }
    # set mappet ((monthly average precipitation)/(monthly average pet)) value to max of 1.25 
    mappet <- ifelse(mappet > 1.25, 1.25, mappet)
    wfac <- (0.2129 + 0.9303 * (mappet) - 0.2413 * (mappet^2))
    # constrain water/moisture effect on decomposition between 0 and 1
    wfac <- ifelse(wfac > 1, 1, ifelse(wfac<0, 0, wfac))
    if(any(IRRIG == 1)){
      wfac[IRRIG == 1] <- wfac_irri
    }
    # calculate temperature effect on decomposition
    tfac <- ((tmax - mtemp)/(tmax - topt))^0.2 * (exp(0.2/2.63 * (1-((tmax - mtemp)/(tmax - topt))^2.63)))
    tfac[is.na(tfac)] <- 0.0001
    # accounting for ground water effect on decomposition (with wfac_gw = 1.5)
    wfac <- wfac*1.5
    # estimate climate decomposition index by combining temperature & water effect on decomposition
    dfac <- mean(wfac)*mean(tfac) 
    dfac <- ifelse(dfac <= 0.0001, 0.0001, dfac)
    #=========================================================================================
    # tillfac - Tillage effect on decomposition (NO-Till is set as reference with tillfac = 1)
    #=========================================================================================
    if(TILL == "FT"){
      tillfac <- tillfac_FT
    }else if(TILL == "RT"){
      tillfac <- tillfac_RT
    }else if(TILL == "NT"){
      # No-Till is the baseline and set to 1
      tillfac <- 1
    }else{
      stop(message("Tillage ", TILL, " not recognized. Must be FT, RT, or NT."))
    }
    #=========================================================================================
    # Calculating f4 and carbon input to the active pool (alpha)
    #=========================================================================================
    f4     <- min(1, max((1.0 - f5 - (0.17 + 0.68 * sand)),0))
    Beta   <- cinput * min(0.999, max((0.85 - 0.018 * (L/N)),0))
    alpha  <- ((Beta * f1) + ((cinput * (1 - L) - Beta) * f2) + (cinput * L * f3 * (f7 + f6 * f8)))/
      (1 - f4 * f7 - f5 * f8 - f4 * f6 * f8)
    #=========================================================================================
    # Calculating Decay Rates (k1, k2, k3, k4, and k5)
    #=========================================================================================
    # decay rate for the structural pool of the litter
    kst <- k10 * dfac
    # decay rate for the metabolic pool of the litter
    kmt <- k20 * dfac * exp((-1 * plig) * L/(1 - (0.85 - 0.018 * (L/N)))) * tillfac
    # decay rate for the active pool
    ka <- k30 * dfac * (0.25 + 0.75 * sand) * tillfac
    # decay rate for the slow pool
    ks <- k40 * dfac * tillfac
    # decay rate for the passive pool
    kp <- k50 * dfac
    #=========================================================================================
    # Calculating Steady State C pools
    #=========================================================================================
    # steady state active pool
    SS_som1   <- alpha/ka
    # steady state slow pool
    SS_som2   <- ((cinput * L) * f3 + (SS_som1 * ka * f4))/ks
    # steady state passive pool
    SS_som3   <- ((SS_som1 * ka * f5) + (SS_som2 * ks * f6))/kp
    #=========================================================================================
    # Adjust decay rates for soil organic matter pools
    # set decay rates to 1 if the rate is > 1 (i.e. residence time < 1 year)
    #=========================================================================================
    # Estimate the carbon stock value for one year time forward
    som1   <- Isom1 + (SS_som1 - Isom1)*ifelse(ka > 1, 1, ka)
    som2   <- Isom2 + (SS_som2 - Isom2)*ifelse(ks > 1, 1, ks)
    som3   <- Isom3 + (SS_som3 - Isom3)*ifelse(kp > 1, 1, kp)
    somsc <- som1 + som2 + som3
    delta.som <- somsc - (Isom1+Isom2+Isom3)
    tempSOM <- tibble("site" = site,"TrtID_Final"=TrtID_Final,"RunIndex"=RunIndex, "year" = year, "cinput_total" = cinput, "dfac" = dfac, 
                          "kA" = ka, "SS_active" = SS_som1, "soc_active" = som1,
                          "kS" = ks, "SS_slow" = SS_som2, "soc_slow" = som2,
                          "kP" = kp, "SS_passive" = SS_som3, "soc_passive" = som3,
                          "soc_total" = somsc, "soc_delta" = delta.som, stringsAsFactors = F)
    SOMstocks <- rbind(SOMstocks, tempSOM)
    # assign new starting value for the next year
    Isom1   <- som1
    Isom2   <- som2
    Isom3   <- som3
  }
  return(SOMstocks)
}