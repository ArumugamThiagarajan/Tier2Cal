## experiment parameters
general_name="test13_Biocluster_ipcc"
date_on="Feb2024"

SIR_sample_size=100
SIR_resample_size=10

## Load all packages
library(tidyr)
library(caret)
library(dplyr)
#library(sensitivity)
library(boot)
library(doParallel)
library(parallel)
library(foreach)
library(purrr)
library(tidyverse)
library(ggplot2)
library(lhs)
library(reshape2)
library(imputeTS)
library(Metrics)
library(here)

## Set main working directory
main.dir <- here::here()
results.dir <- gsub("aafc-soc-cal","aafc-soc-cal-results",main.dir)


## Load data
site_data <- read.csv(paste0(main.dir,"//data//input data//site data//LTE_Master_beta_nov_15.csv"))
#dsm_data<- read.csv(paste0(main.dir,"//data//input data//site data//lte_soc_30_cm_dsm_expid_may_2.csv"))
climate.dir=paste0(main.dir,"//data//input data//climate data//W9param_TablesCleaned")
parameter_bounds <- read.csv(paste0(main.dir,"//data//input data//parameter data//ipcct2_parameters_gsa.csv"), stringsAsFactors = FALSE)
cinput.df<-read.csv(paste0(main.dir,"//data//prepared data//yield_and_cinput_prepared.csv"))
polyids_all <- read.csv(paste0(main.dir,"//data//input data//c_input//exp_location_SLC_linkage.csv")) %>%
  select(Exp_ID,POLYID)


## Load functions
#source(paste0(main.dir,"//scripts//Recent//Top30cm_SOC_estimator.R")) # Tool to get standardized SOC values for top 30 cm of soil
#source(paste0(main.dir,"//scripts//Recent//data preparation ipcct2.R")) # Data preperation function
source(paste0(main.dir,"//scripts//Recent//ipcct2_SOCTier2Model.r")) # IPCC tier 2 steady state model
source(paste0(main.dir,"//scripts//Recent//Spinnup IPCCT2.R")) # Model initiation for IPCCt2
#source(paste0(main.dir,"//scripts//Recent//sensitivity analysis IPCCT2.R")) # GSA function for IPCCt2
source(paste0(main.dir,"//scripts//Recent//SIR IPCCT2.R")) # SIR function to obtain posterior distribution
#source(paste0(main.dir,"//scripts/Recent/create run map.R")) # Treatment-stocks map to run the experiment
source(paste0(main.dir,"//scripts/Recent/predictor ippct2.R")) # Final preparation step for the ipcct2 model
source(paste0(main.dir,"//scripts/Recent/Summary fit statistics.r")) # Get summary fit statistics from observed vs predicted values
source(paste0(main.dir,"//scripts/Recent/run.ipcct2.par.R")) # Generic function to run the IPCCt2 model by runindex rapidly

gsa_sample_size=1000

Experiment_name=paste0(general_name,"_",date_on,"_GSA ",gsa_sample_size,"_SIR ",SIR_sample_size,"To",SIR_resample_size)

folders_needed=c(paste0(main.dir,"//data//prepared data//",Experiment_name),
                 paste0(main.dir,"//data//prepared data//",Experiment_name,"//IPCCT2"),
                 paste0(main.dir,"//results//",Experiment_name),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//GSA"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//SIR"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//SIR//interim"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//validation//"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//prior//"),
                 paste0(main.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior//"))

for (i in folders_needed){
  ifelse(dir.exists(i),FALSE, dir.create(i))
}


## Create folder structure

folders_needed=c(paste0(main.dir,"//data//prepared data"),
                 paste0(main.dir,"//data//prepared data//IPCCT2"),
                 paste0(main.dir,"//results"),
                 paste0(main.dir,"//results//IPCCT2"),
                 paste0(main.dir,"//results//IPCCT2//spin up"),
                 results.dir,
                 paste0(results.dir,"//results//"),
                 paste0(results.dir,"//results//",Experiment_name),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//GSA"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//interim"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//prior//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior all calibration//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior all validation//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//prior all validation//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//prior values validation//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//prior values calibration//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior values validation//"),
                 paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//validation//interim//posterior values calibration//"))

for (i in folders_needed){
  ifelse(dir.exists(i),FALSE, dir.create(i))
}

### Prepared data load

stocks.df <-read.csv(paste0(main.dir,"//data//prepared data//IPCCT2//modelled_stocks_top_30_cm_nov20_2023.csv"))
site_data_ipcct2=read.csv(paste0(main.dir,"//data//prepared data//IPCCT2//sites_prepared_ipcct2.csv"))
climate_data_ipcct2 <- read.csv(paste0(main.dir,"//data//prepared data//IPCCT2//climate_prepared_ipcct2.csv"))

RunMap<- read.csv(paste0(main.dir,
                         "//data//prepared data//SOC-cal_RunMap_with_initSOC_Nov2023_NP.csv"))%>%
  filter(valid==1)

# Default parameters for IPCCt2
params.default=data.frame(tillfac_FT = 3.036,
                          tillfac_RT = 2.075,
                          wfac_irri  = 0.775,
                          k10        = 18.5,
                          k20        = 4.9,
                          k30        = 7.4,
                          k40        = 0.209,
                          k50        = 0.00689,
                          f1         = 0.378,
                          f2         = 0.368,
                          f3         = 0.455,
                          f5         = 0.0855,
                          f6         = 0.0504,
                          f7         = 0.42,
                          f8         = 0.45,
                          tmax       = 45,
                          topt       = 33.69,
                          plig       = 3)



### Split calibration and evaluation datasets
# Balanced split for different drivers
calibration_sites <- RunMap%>% filter(group!=1)
evaluation_sites <- RunMap%>% filter(group==1)

write.csv(calibration_sites,file=paste0(main.dir,
                                        "//data//prepared data//",
                                        Experiment_name,
                                        "//RunMap_cal.csv"),row.names=F)
write.csv(evaluation_sites,file=paste0(main.dir,
                                       "//data//prepared data//",
                                       Experiment_name,
                                       "//RunMap_val.csv"),row.names=F)


### Create the final site predictor datasets for IPCCt2
calibration_site_predictor<-predictor.ipcct2(df=calibration_sites,
                                             prepared.data=site_data_ipcct2)
evaluation_site_predictor<-predictor.ipcct2(df=evaluation_sites,
                                            prepared.data=site_data_ipcct2)

# Save to file
write.csv(calibration_site_predictor,file=paste0(main.dir,
                                                 "//data//prepared data//",
                                                 Experiment_name,
                                                 "//IPCCT2//site_data_calibration.csv"),row.names=F)
write.csv(evaluation_site_predictor,file=paste0(main.dir,
                                                "//data//prepared data//",
                                                Experiment_name,
                                                "//IPCCT2//site_data_evaluation.csv"),row.names=F)


## Sensitivity analysis output
sensitive_params<- data.frame(params=c("f1","f5","f6","k50","tmax","topt","tillfac_RT","tillfac_FT"))

new_params_bounds <- parameter_bounds %>% 
  mutate(cal=ifelse(Parameter %in% sensitive_params$params,1,0)) %>%
  mutate(lower=ifelse(cal==0,value,lower),
         upper=ifelse(cal==0,value,upper))


# SIR calibration
## Make the LHS sample

#=================================================================================
# read prior distribution
# (Required columns: Parameter, value, lower, upper)
paramBounds <- parameter_bounds
#=================================================================================
# names of parameters that are allowed to vary
varSI       <- sensitive_params$params
nParams     <- length(varSI)
#=================================================================================
# LHS sampling for SIR
#=================================================================================
# sample size (1000000 samples were used in Gurung et al., 2020)
n <- SIR_sample_size

set.seed(123)
X1 <- randomLHS(n = n, k = nParams)
# transform standard uniform LHS to prior distribution 
Y1 <- matrix(0, nrow = nrow(X1), ncol = nParams)
for(i in 1:nParams){
  pos <- which(paramBounds$Parameter == varSI[i])
  lower <- paramBounds[pos, "lower"]
  upper <- paramBounds[pos, "upper"]
  Y1[, i] <- qunif(X1[, i], min = lower, max = upper)
}
X <- as.data.frame(Y1)
names(X) <- varSI
X <- cbind("SampleID" = 1:nrow(X), X)

names_X <- names(X)
for (i in 1:nrow(parameter_bounds)){
  param_on <- parameter_bounds[i,1]
  param_def <- parameter_bounds[i,2]
  if (param_on %in% names_X){} else 
  {X[[param_on]]<-param_def}
}

write.csv(X,paste0(main.dir,
                   "//data//prepared data//",
                   Experiment_name,
                   "//IPCCT2//prior_parameter_distribution_LHS.csv"),row.names=F)


## Perform SIR calibration on IPCCt2 model parameters
calibration_data=calibration_site_predictor

run.IPCCT2.par(main.dir=main.dir,
               Experiment_name = Experiment_name,
               results.dir = results.dir,
               out.subfolder="SIR//interim//",
               parameter.df=X,
               stocks_data=stocks.df,
               calibration.df=calibration_data,
               climate_data_ipcct2 = climate_data_ipcct2,
               site_data_ipcct2=site_data_ipcct2,
               RunMap=RunMap,
               out_type = "fit",
               likelihood="delta")

loglike_out <- list.files(path = paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//interim//"),full.names = T)

readFunction <- function (x) {
  return(tryCatch(read.csv(x), error=function(e) NULL))
}

tables <- lapply(loglike_out, readFunction)

combined.df <- do.call(rbind , tables)

Lkhood1 <- combined.df %>%
  mutate(loglike = ifelse(loglike == -Inf, NA, loglike)) %>%
  group_by(SampleID) %>%
  summarise(loglike = mean(loglike, na.rm=T),
            RMSE = mean(RMSE, na.rm=T),
            bias = mean(bias,na.rm=T)) %>%
  mutate(weights = exp(loglike)/sum(exp(loglike)))

write.csv(Lkhood1 ,
          paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//prior_calibration_fit.csv"),row.names=F)

Lkhood1<- read.csv(paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//prior_calibration_fit.csv"))
sampIndx <- sample(1:nrow(Lkhood1), 
                   size = SIR_resample_size,
                   replace = FALSE,
                   prob = Lkhood1$weights)

X<-read.csv(paste0(main.dir, "//data//prepared data//",
                   Experiment_name,
                   "//IPCCT2//prior_parameter_distribution_LHS.csv"))

PostTheta <- as.data.frame(X[sampIndx,])

# X <- combined.df %>% select(-V1,-RunIndex,-loglike,-RMSE) %>% group_by(SampleID) %>% summarise_all(mean) # if prior was not saved

PostTheta <- X %>%
  filter(SampleID %in% sampIndx)


write.csv(PostTheta ,
          paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//SIR_posterior_distribution.csv"),row.names=F)


# Find maximum a posteriori
get_max_density <- function(x) {
  dens <- x %>%
    density
  return(dens$x[which.max(dens$y)])
}

sir_table <- PostTheta %>%
  reshape2::melt(id.vars="SampleID") %>%
  group_by(variable)%>%
  summarise(`2.5%` = quantile(value, probs = 0.025),
            `25%` = quantile(value, probs = 0.25),
            median = quantile(value, probs = 0.50),
            `75%` = quantile(value, probs = 0.75),
            `97.5%` = quantile(value, probs = 0.975),
            map = get_max_density(value),
            mean = mean(value)) %>%
  rename(Parameter = variable) %>% 
  merge(parameter_bounds %>% select(Parameter,default=value,lower,upper),by="Parameter")

write.csv(sir_table ,paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//SIR_summary.csv"),row.names=F)

## Graphing:

prior <- X %>%
  reshape2::melt(id.vars="SampleID") %>%
  mutate(distribution="prior")

posterior <- PostTheta %>%
  reshape2::melt(id.vars="SampleID") %>%
  mutate(distribution="posterior")

df.plot <- rbind(prior,posterior) %>%
  filter(variable %in% sensitive_params$params)

sir_graph <- ggplot() +
  geom_density(data = df.plot, aes(value,fill=distribution), col = NA, alpha = 0.5) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_classic()

ggsave(plot = sir_graph, paste0(results.dir,"//results//",Experiment_name,"//IPCCT2//SIR//SIR_distributions.png"))