#Plankton compilation - data processing

library(dplyr)
library(tidyr)
source("scripts/chaos_localstability_methods.R")

#files from Rogers et al. 2023
site_meta=read.csv("data/plankton/original/sitemetadata.csv")
plankton_meta=read.csv("data/plankton/original/seriesmetadata.csv")
plankton_ts=read.csv("data/plankton/original/alldata_filtered.csv")
plankton_results=read.csv("data/plankton/original/results_global_updated.csv")

plankton_tsn=plankton_ts %>% group_by(Site,Name,Level) %>% nest(.key="data_rescale") %>%
  mutate(data_rescale=map(data_rescale, as.data.frame))
plankton=left_join(plankton_tsn, plankton_meta, by = c("Site", "Name", "Level"))
plankton=left_join(plankton_results %>% select(Site:theta,JLEmin,JLEmean,JLEsign), plankton, by = c("Site", "Name", "Level"))
plankton=left_join(plankton, site_meta %>% select(Site,Country,Type,Redistributable,Lat,Lon))

plankton=plankton %>% filter(Level=="SP", longestrun>=40, Site!="Lake Kasumigaura")

plankton_full=unnest(plankton, data_rescale)

#time series
plankton_series=select(plankton_full, colnames(plankton_ts))
plankton_series_public=filter(plankton_full, Redistributable=="yes") %>%
  select(colnames(plankton_ts))
unique(plankton_series$Site)
unique(plankton_series_public$Site)

#metadata
plankton_seriesmeta=select(plankton, colnames(plankton_meta), Site,Country,Type,Redistributable,Lat,Lon)

#results
plankton_seriesresults=plankton %>%
  mutate(JLEmean_mo=JLEmean, taumax=12) %>%
  select(Site:R2gr,modelform:theta,JLEmin,JLEmean,JLEmean_mo,JLEsign,taumax)

#add body size information
#both these tables from Rogers et al. 2022
# gpdd_d=read.csv("data/GPDD/original/gpdd_ts_metadata.csv")
# plank=read.csv("data/plankton/original/lakes_ts_metadata.csv", stringsAsFactors = F)
# planktraits=gpdd_d %>% select(Name=TaxonName, Mass_g, MinAge_mo) %>%
#   rbind(plank %>% select(Name, Mass_g, MinAge_mo))
# plankton_lh= plankton_seriesmeta %>% select(Name, Site) %>%
#   mutate(Binomial2=gsub("[._]"," ", Name)) %>%
#   left_join(planktraits, by = c("Binomial2"="Name")) %>%
#   left_join(planktraits, by = c("Name"="Name")) %>%
#   arrange(Binomial2)
# write.csv(plankton_lh, "data/plankton/plankton_lifehistory.csv", row.names = F)
# manually edit. these were also filtered to only those that met the post-analysis
# criteria for inclusion (e.g. R2>0.15) to reduce the number of species we had to
# find information for
plankton_lh=read.csv("data/plankton/plankton_lifehistory_manual.csv", stringsAsFactors = F)
plankton_seriesmeta=left_join(plankton_seriesmeta, plankton_lh %>% select(Name,Site,Mass_g,MinAge_mo)) %>%
  mutate(SamplingInterval="monthly", mo_timestep=1)

#get monotonic trend (not previously computed)
monotonic_eval=function(data) {
  data$decYear=data$Year+(data$Month-1)/12
  cor(data$Year, data$Abundance, use="p", method="spearman")
}
plankton_seriesmeta$monotonic=map_dbl(plankton$data_rescale, monotonic_eval)

#add memory stats to results (not previously computed)
load("data/plankton/original/ts_results_updated3.Rdata") #contains jacobians
ts_results2$memorystats=map(ts_results2$jacobians, memorystats)
memstats=unnest(ts_results2 %>% select(Site, Name, Level, memorystats), cols = memorystats)

plankton_seriesresults=plankton_seriesresults %>%
  left_join(memstats %>% select(Site, Name, Level, meanlag))

#write new files
write.csv(plankton_series,"data/plankton/plankton_series.csv", row.names = F)
write.csv(plankton_series_public,"data/plankton/plankton_series_public.csv", row.names = F)
write.csv(plankton_seriesmeta,"data/plankton/plankton_seriesmetadata.csv", row.names = F)
write.csv(plankton_seriesresults,"results/plankton_LEresults.csv", row.names = F)
