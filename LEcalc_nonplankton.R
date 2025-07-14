# Applies Jacobian LE chaos detection method to nonplankton series
# GPDD, LP, InsectChange

# LEs were already computed for GPDD in Rogers et al. 2022, but redoing with
# maxtau=12 for monthly series for consistency with other analyses. For
# annual series, maxtau=6, results will be identical. The results do not
# change substantially.

library(rEDM)
library(dplyr)
library(tidyr)
library(purrr)

source("scripts/chaos_localstability_methods.R")

ts1_d=read.csv("data/GPDD/GPDD_series.csv") %>% group_by(Id) %>% nest(.key="data_rescale")
ts1_meta=read.csv("data/GPDD/GPDD_seriesmetadata.csv") %>%
  mutate(Source="GPDD") %>% select(Source, Id, Loc, Name, SamplingInterval, mo_timestep)

ts2_d=read.csv("data/LP/LP_series.csv") %>% group_by(Id) %>% nest(.key="data_rescale")
ts2_meta=read.csv("data/LP/LP_seriesmetadata.csv") %>%
  mutate(Source="LP") %>% select(Source, Id, Loc, Name, SamplingInterval, mo_timestep)

ts3_d=read.csv("./data/insects/insects_series.csv") %>% group_by(Id) %>% nest(.key="data_rescale")
ts3_meta=read.csv("data/insects/insects_seriesmetadata.csv") %>%
  mutate(Source="InsectChange") %>% select(Source, Id, Loc, Name, SamplingInterval, mo_timestep)

#stack
ts_d=rbind(ts1_d,ts2_d,ts3_d)
ts_meta=rbind(ts1_meta,ts2_meta,ts3_meta)

ts_d=ts_d %>%
  mutate(data_rescale=map(data_rescale, as.data.frame)) %>%
  left_join(ts_meta) %>%
  mutate(taumax=ifelse(SamplingInterval=="monthly", 12, 6))

ts_results=select(ts_d, Id, Loc, Name)

#for updating
#dnestfilt=ts_d

#### Jacobian LE Method ####

#fit models
#models 1 and 2 are dropped because the results of 1&3 and 2&5 are identical
ts_results$modelresults3=pmap(list(data=ts_d$data_rescale, taumax=ts_d$taumax), smap_model_options, y="AbundanceScale", model=3, seaspred=F) #fd-ut
ts_results$modelresults4=pmap(list(data=ts_d$data_rescale, taumax=ts_d$taumax), smap_model_options, y="AbundanceScale", model=4, seaspred=F) #gr-ut
ts_results$modelresults5=pmap(list(data=ts_d$data_rescale, taumax=ts_d$taumax), smap_model_options, y="AbundanceScale", model=5, seaspred=F) #gr-log

#pull out R2 values for each model
ts_d$R3m=map_dbl(ts_results$modelresults3, ~.x$modelstats$R2model)
ts_d$R4m=map_dbl(ts_results$modelresults4, ~.x$modelstats$R2model)
ts_d$R5m=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2model)
ts_d$R3a=map_dbl(ts_results$modelresults3, ~.x$modelstats$R2abund)
ts_d$R4a=map_dbl(ts_results$modelresults4, ~.x$modelstats$R2abund)
ts_d$R5a=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2abund)

#get best model of the 3 forms
ts_d$bestmodel=select(ts_d %>% ungroup(),R3a,R4a,R5a) %>% apply(1,which.max)
ts_d$bestR2=select(ts_d %>% ungroup(),R3a,R4a,R5a) %>% apply(1,max) #R2 for abundance
ts_d$bestR2m=select(ts_d %>% ungroup(),R3m,R4m,R5m) %>% apply(1,max) #R2 for growth rate
ts_results$modelresultsbest=cbind(select(ts_results %>% ungroup(), modelresults3, modelresults4, modelresults5),ts_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["ts_d$bestmodel"]); x[m][[1]]})
#get E, tau, theta values from best model
ts_d$E=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$E)
ts_d$tau=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$tau)
ts_d$theta=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$theta)
#get best model form
ts_d$modelform=map_chr(ts_results$modelresultsbest, ~.x$form)

#get jacobian matrices
ts_results$jacobians=map(ts_results$modelresultsbest, getJacobians)

#get global LE
ts_results$GlobalLE=map2(ts_results$modelresultsbest, ts_results$jacobians, LEshift)
ts_d$JLEmin=map_dbl(ts_results$GlobalLE, ~.x$minci) #LE lower confidence bound (*this is the LE estimate to use!*)
ts_d$JLEmean=map_dbl(ts_results$GlobalLE, ~.x$minmean) #LE mean (*this is the LE estimate to use!*)
ts_d$JLEsign=ifelse(ts_d$JLEmin>0.01, "chaotic", "not chaotic")
ts_d$JLEmean_mo=ts_d$JLEmean/ts_d$mo_timestep

#get memory stats
ts_results$memorystats=map(ts_results$jacobians, memorystats)

#### Export Results ####

#save results
save(ts_d, ts_results, file = "results/nonplankton_LE_07-09-25.Rdata")

#global LE and memory stats
memstats=unnest(ts_results %>% select(Id, memorystats), cols = memorystats)
exportres1=ts_d %>% ungroup() %>%
  select(Id, Loc, Name, R2abund=bestR2, R2gr=bestR2m, modelform, E, tau, theta,
         JLEmin, JLEmean, JLEmean_mo, JLEsign, taumax) %>%
  left_join(memstats %>% select(Id, meanlag))
write.csv(exportres1, "results/nonplankton_LEresults.csv", row.names = F)
