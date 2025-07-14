# Analysis of nonplankton series (fit of alternative models)
# GPDD, LP, InsectChange

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rEDM)
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

#models ####

#loglinear
ts_d$loglinear=map(ts_d$data_rescale, fit_loglinear, var="AbundanceScale", exrad=0)
ts_d$loglinearR2=map_dbl(ts_d$loglinear, function(d) GPEDM::getR2(d$obs, d$pred))

#ricker
ts_d$ricker=map(ts_d$data_rescale, fit_ricker, var="AbundanceScale", exrad=0)
ts_d$rickerR2=map_dbl(ts_d$ricker, function(d) GPEDM::getR2(d$obs, d$pred))

#1d smap (memoryless)
ts_d$smap1d=pmap(list(data=ts_d$data_rescale, taumax=ts_d$taumax), best1dsmap, y="AbundanceScale")
ts_d$smap1dR2=map_dbl(ts_d$smap1d, ~.x$modelstats$R2abund)

#best of 1d models
ts_d$best1dR2=select(ts_d %>% ungroup(),loglinearR2,rickerR2,smap1dR2) %>% apply(1,max)
ts_d$best1dmodel=select(ts_d %>% ungroup(),loglinearR2,rickerR2,smap1dR2) %>% apply(1,which.max)
ts_d$best1dmodel=ifelse(ts_d$best1dmodel==1,"loglinear", ifelse(ts_d$best1dmodel==2, "ricker", "1dsmap"))
table(ts_d$best1dmodel)

#best AR (linear)
ts_d$ARmod=pmap(list(data=ts_d$data_rescale, taumax=ts_d$taumax), bestAR, y="AbundanceScale")
ts_d$ARmodR2=map_dbl(ts_d$ARmod, ~.x$modelstats$R2abund)

#save output
save(ts_d, file = "results/nonplankton 07-09-25.Rdata")

ts_out=ts_d %>% select(Id, Loc, Name,
                       smap1dR2, loglinearR2, rickerR2,
                       best1dR2, best1dmodel, ARmodR2)
write.csv(ts_out, "results/nonplankton_results.csv", row.names = F)
