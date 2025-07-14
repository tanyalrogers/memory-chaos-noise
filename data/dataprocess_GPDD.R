# GPDD - data processing

# Subsetting to series with at least 40 contiguous nonmissing timepoints
# Also excluding Port Erin Bay plankton, since included in the plankton dataset

library(dplyr)
library(tidyr)
library(purrr)

#files from Rogers et al 2022
ts_d=read.csv("data/GPDD/original/gpdd_timeseries.csv")
ts_meta=read.csv("data/GPDD/original/gpdd_ts_metadata.csv")

ts_d=ts_d %>%
  #rename cols for consistency
  mutate(Id=paste('GPDD',MainID, sep = "_")) %>%
  select(Id, SeriesStep, SampleYear, Abundance=PopulationUntransformed,
         AbundanceScale=PopRescale, AbundanceLogScale=PopRescale_log) %>%
  group_by(Id) %>% nest(.key="data_rescale") %>%
  mutate(data_rescale=map(data_rescale, as.data.frame))
ts_meta=ts_meta %>%
  mutate(Id=paste("GPDD",MainID, sep="_"), Loc=substr(ExactName,1,40), Name=CommonName) %>%
  relocate(Id,Loc,Name) %>%
  mutate(monotonic=sqrt(monotonicR2))

#subset time series
ts_meta$longestrun=map_dbl(ts_d$data_rescale, function(d) {
  len=rle(!is.na(d$AbundanceScale))
  len2=len$lengths[len$values==TRUE]
  max(len2)
})
ts_meta=subset(ts_meta, longestrun>=40 & LocationID!=1024) #remove Port Erin Bay plankton data

#subset cols
ts_meta=ts_meta %>%
  select(Id,Loc,Name,MainID:SamplingInterval,mo_timestep,ndatapoints,longestrun,propmissing,propzeros,monotonic,
         Mass_g,MinAge_mo,TaxonomicClass3,ExactName,Country,LongDD,LatDD,data_rescale_case)

ts_d=ts_d %>%
  filter(Id %in% ts_meta$Id)
ts_out=unnest(ts_d,data_rescale) %>% ungroup()

#write new files
write.csv(ts_out,"data/GPDD/GPDD_series.csv", row.names = F)
write.csv(ts_meta,"data/GPDD/GPDD_seriesmetadata.csv", row.names = F)
