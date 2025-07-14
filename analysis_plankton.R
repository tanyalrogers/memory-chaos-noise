# Analysis of Plankton series (fit of alternative models)
# all series are monthly

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rEDM)
source("scripts/chaos_localstability_methods.R")

plankton=read.csv("data/plankton/plankton_seriesmetadata.csv")
plankton_ts=read.csv("data/plankton/plankton_series.csv")

plankton_ts=plankton_ts %>% group_by(Site,Name,Level) %>% nest(.key="data_rescale") %>%
  mutate(data_rescale=map(data_rescale, as.data.frame))
plankton=left_join(plankton, plankton_ts, by = c("Site", "Name", "Level"))

#models ####

#loglinear
plankton$loglinear=map(plankton$data_rescale, fit_loglinear, var="AbundanceScale", exrad=0)
plankton$loglinearR2=map_dbl(plankton$loglinear, function(d) GPEDM::getR2(d$obs, d$pred))

#ricker
plankton$ricker=map(plankton$data_rescale, fit_ricker, var="AbundanceScale", exrad=0)
plankton$rickerR2=map_dbl(plankton$ricker, function(d) GPEDM::getR2(d$obs, d$pred))

#1d smap (memoryless)
plankton$smap1d=map(plankton$data_rescale, best1dsmap, y="AbundanceScale")
plankton$smap1dR2=map_dbl(plankton$smap1d, ~.x$modelstats$R2abund)

#best of 1d models
plankton$best1dR2=select(plankton %>% ungroup(),loglinearR2,rickerR2,smap1dR2) %>% apply(1,max)
plankton$best1dmodel=select(plankton %>% ungroup(),loglinearR2,rickerR2,smap1dR2) %>% apply(1,which.max)
plankton$best1dmodel=ifelse(plankton$best1dmodel==1,"loglinear", ifelse(plankton$best1dmodel==2, "ricker", "1dsmap"))
table(plankton$best1dmodel)

#best AR (linear)
plankton$ARmod=map(plankton$data_rescale, bestAR, y="AbundanceScale")
plankton$ARmodR2=map_dbl(plankton$ARmod, ~.x$modelstats$R2abund)

#save output
save(plankton, file = "results/plankton 07-09-25.Rdata")

plankton_out=plankton %>% select(Site, Name, Level,
                                 smap1dR2, loglinearR2, rickerR2,
                                 best1dR2, best1dmodel, ARmodR2)
write.csv(plankton_out, "results/plankton_results.csv", row.names = F)
