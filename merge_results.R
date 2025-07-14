# Merge results and metadata tables from all datasets
# Filter series (no exponential growth, R2>0.15, no strong monotonic trends)
# Compute some additional stats

library(dplyr)
library(tidyr)

# load data ####

#full model, LE, and memory results
nonplankton_full=read.csv("results/nonplankton_LEresults.csv", stringsAsFactors = F)
plankton_full=read.csv("results/plankton_LEresults.csv", stringsAsFactors = F)

#alternative model results
nonplankton_alt=read.csv("results/nonplankton_results.csv", stringsAsFactors = F)
plankton_alt=read.csv("results/plankton_results.csv", stringsAsFactors = F)

#metadata
gpdd_meta=read.csv("data/GPDD/GPDD_seriesmetadata.csv", stringsAsFactors = F)
lp_meta=read.csv("data/LP/LP_seriesmetadata.csv", stringsAsFactors = F)
insects_meta=read.csv("data/insects/insects_seriesmetadata.csv", stringsAsFactors = F)
plankton_meta=read.csv("data/plankton/plankton_seriesmetadata.csv", stringsAsFactors = F)

#standardize metadata ####

gpdd_meta=gpdd_meta %>% select(Id,Loc,Name,SamplingInterval,mo_timestep,ndatapoints,longestrun,monotonic,Mass_g,TaxonomicClass3) %>%
  mutate(Group="Non-plankton species", Source="GPDD",
         TaxonomicGroup=ifelse(TaxonomicClass3 %in% c("Birds","Mammals"),"Birds/Mammals",
                               ifelse(TaxonomicClass3=="Bony fishes","Fishes", TaxonomicClass3)))
lp_meta=lp_meta %>% select(Id,Loc,Name,SamplingInterval,mo_timestep,ndatapoints,longestrun,monotonic,Mass_g,TaxonomicClass3=TaxonomicGroup) %>%
  mutate(Group="Non-plankton species", Source="LP",
         TaxonomicGroup=ifelse(TaxonomicClass3 %in% c("Birds","Mammals"),"Birds/Mammals",TaxonomicClass3))
insects_meta=insects_meta %>% select(Id,Loc,Name,SamplingInterval,mo_timestep,longestrun,ndatapoints,monotonic,Mass_g) %>%
  mutate(Group="Non-plankton species", Source="InsectChange", TaxonomicClass3="Insects", TaxonomicGroup="Insects")
plankton_meta=plankton_meta %>% select(Site,Name,Level,SamplingInterval,mo_timestep,ndatapoints,longestrun,monotonic,Mass_g) %>%
  mutate(Group="Plankton species", Source="PlanktonCompilation", TaxonomicClass3="Plankton", TaxonomicGroup="Plankton")

#merge tables ####

nonplankton_meta=rbind(gpdd_meta, lp_meta, insects_meta)
nonplankton_join=nonplankton_meta %>% left_join(nonplankton_full) %>% left_join(nonplankton_alt)
plankton_join=plankton_meta %>% left_join(plankton_full) %>% left_join(plankton_alt) %>%
  mutate(Id=paste("P",1:nrow(plankton_meta),sep="_")) %>%
  rename(Loc=Site) %>% relocate(Id,Loc,Name) %>% select(-Level)
allresults=rbind(nonplankton_join, plankton_join) %>% rename(smapoptR2=R2abund)

allresults$Group=factor(allresults$Group, levels = unique(allresults$Group))
allresults$SamplingInterval2=ifelse(allresults$SamplingInterval=="monthly","monthly","annual")
table(allresults$Group)
table(allresults$TaxonomicGroup)
table(allresults$TaxonomicClass3)

#filter ####

#remove series with low predictability, strong monotonic trends, or exponential growth (linear and positive LE)
allresults$expgrowth=ifelse(allresults$theta==0 & allresults$JLEsign=="chaotic", "yes", "no")
allresults_filt00=allresults %>% filter(smapoptR2>0.15)
allresults_filt0=allresults_filt00 %>% filter(abs(monotonic)<0.8 & expgrowth=="no")

#removed due to low predictability
allresults %>% filter(smapoptR2<=0.15) %>% nrow()
#removed due to monotonic trend or exp growth
allresults_filt00 %>% filter(abs(monotonic)>=0.8 | expgrowth=="yes") %>% nrow()
#series with exp growth
allresults_filt00 %>% filter(expgrowth=="yes") %>% nrow()
allresults_filt00 %>% filter(abs(monotonic)>=0.8) %>% nrow()
allresults_filt00 %>% filter(abs(monotonic)>=0.8 & expgrowth=="yes") %>% nrow()
allresults_filt00 %>% filter(abs(monotonic)>=0.8 & smapoptR2<=0.15) %>% nrow()

#removed=allresults_filt00 %>% filter(abs(monotonic)>=0.8 | expgrowth=="yes")

#final filtered series
allresults_filt=allresults_filt0

#compute more values ####

allresults_filt$diffbest1d=allresults_filt$smapoptR2-allresults_filt$smap1dR2
allresults_filt$diffAR=allresults_filt$smapoptR2-allresults_filt$ARmodR2
allresults_filt$JLEsign2=ifelse(abs(allresults_filt$JLEmean)<0.02, "neutral", ifelse(allresults_filt$JLEsign=="chaotic", "unstable", "stable"))
allresults_filt$JLEsign2_mo=ifelse(abs(allresults_filt$JLEmean_mo)<0.02, "neutral", ifelse(allresults_filt$JLEsign=="chaotic", "unstable", "stable"))
allresults_filt$JLEsign2=factor(allresults_filt$JLEsign2, levels = c("stable", "neutral", "unstable"))
allresults_filt$JLEsign2_mo=factor(allresults_filt$JLEsign2_mo, levels = c("stable", "neutral", "unstable"))

#standardize memory stats
allresults_filt = allresults_filt %>%
  mutate(meanlag_st=meanlag*tau*mo_timestep)

#output combined data ####

write.csv(allresults_filt, "results/combined_results.csv", row.names = F)
