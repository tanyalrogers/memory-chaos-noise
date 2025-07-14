#Living Planet time series - data processing

library(dplyr)
library(tidyr)
library(purrr)

#downloaded data file (converted to UTF-8 encoding)
dtimeseries=read.csv("C:/Users/trogers/Documents/Nonstationarity/LPD2022_public/LPD2022_public2.csv",
                     na.strings = "NULL")

dmeta=select(dtimeseries, 1:33)
dnest0=dtimeseries %>% pivot_longer(34:104, names_to = "Year", values_to = "Abundance") %>%
  mutate(Year=as.numeric(sub("X","",Year)))
dnest0=dnest0 %>% select(ID, Year, Abundance) %>% group_by(ID) %>% nest()

dnest0$ndatapoints=map_dbl(dnest0$data, function(d) {
  length(which(!is.na(d$Abundance)))})
#at least 40 data points
dnest=filter(dnest0, ndatapoints>40)

#trim off and remove long periods with no dynamics
tsprocess=function(data) {
  # AbundancePos=ifelse(data$Abundance<=0,NA,data$Abundance)
  # for(i in 1:(length(AbundancePos)-12)) {
  #   if(all(is.na(AbundancePos[i:(i+12)]))) {
  #     data$Abundance[i:(i+12)]<-NA
  #   }
  # }
  # #data$AbundanceScale=data$AbundancePos/sd(data$AbundancePos, na.rm = T)
  # #data$AbundanceLogScale=log(data$AbundanceScale)
  data=data[(first(which(!is.na(data$Abundance)))):(last(which(!is.na(data$Abundance)))),]
  return(as.data.frame(data))
}
dnest$data2=map(dnest$data,tsprocess)

#figure out what to do with zeros
# dnest$mina=map_dbl(dnest$data2, ~min(.x$Abundance[which(.x$Abundance>0)], na.rm=T))
# dnest$minb=map_dbl(dnest$data2, ~min(.x$Abundance, na.rm=T))
# dnest$minc=round(dnest$mina-dnest$minb,2)
# dnest$maxa=map_dbl(dnest$data2, ~max(.x$Abundance, na.rm=T))
# dnest=dnest %>% arrange(desc(minc))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
rescale=function(data, diag=F) {
  mind=min(data$Abundance, na.rm=T)
  minnz=min(data$Abundance[data$Abundance>0], na.rm=T)
  if(mind>0) {
    data$AbundanceScale=data$Abundance/sd(data$Abundance, na.rm=T)
    case=1 #no zeros, nothing added
  } else {
    if(all(is.wholenumber(data$Abundance), na.rm=T) & minnz<100) {
      data$AbundanceScale=(data$Abundance+1)/sd(data$Abundance, na.rm=T)
      case=2 #zeros, integers, plus 1, unless minnz is >100
    } else {
      data$AbundanceScale=(data$Abundance+minnz)/sd(data$Abundance, na.rm=T)
      case=3 #zeros, non-integers, plus min non-zero value
    }
  }
  data$AbundanceLogScale=log(data$AbundanceScale)
  if(diag) return(case)
  else return(as.data.frame(data))
}
dnest$data_rescale=map(dnest$data2, rescale)
dnest$data_rescale_case=map_dbl(dnest$data2, rescale, diag=T)

#filtering metrics
dnest$tslength=map_int(dnest$data_rescale,nrow)

dnest$longestrun=map_int(dnest$data_rescale, function(d) {
  len=rle(!is.na(d$Abundance))
  max(len$lengths[len$values==TRUE])
})
dnest$propmissing=map_dbl(dnest$data_rescale, function(d) {
  ts=d$Abundance #[(first(which(!is.na(d$AbundancePos)))):(last(which(!is.na(d$AbundancePos))))]
  length(which(is.na(ts)))/length(ts)
})
dnest$propzeros=map_dbl(dnest$data_rescale, function(d) {
  ts=d$Abundance #[(first(which(!is.na(d$AbundancePos)))):(last(which(!is.na(d$AbundancePos))))]
  length(which(ts==0))/length(ts)
})
#number of unique values
dnest$uniquevals=map_dbl(dnest$data_rescale, function(d) {
  n_distinct(d$Abundance, na.rm = T)
})

#filter
dnestfilt=dnest %>%
  filter(longestrun>=40) %>%
  filter(propzeros<0.6) %>%
  filter(propmissing<0.22) %>%
  filter(uniquevals>=5)

#exclude duplicates in gpdd (see section below for how these were identified)
exclude=c(27571,1047,6053,6098, 22059, 22062,22063,22073,22076)
dnestfilt=dnestfilt %>%
  filter(!(ID %in% exclude))
dmeta2=dmeta %>% filter(ID %in% dnestfilt$ID)

#remove 'data' that are actually model output

#population models (mostly fish, some marine mammals)
removed=filter(dmeta2,grepl(c("Integrated Analysis|Assessed|Stock assessment|XSA|Stock Synthesis|catch at age|catch-at-age|Population model|AD-CAM|Biomass dynamics model|surplus production model|age-based assessment|age-structured assessment|spawning stock"), dmeta2$Method, ignore.case = T))
unique(removed$Method)
#breeding bird survey (model of annual deviations imposed on growth trend)
#look mostly like exp growth/decay because this behavior is built in to model. not useful for our purposes
removed2=filter(dmeta2,grepl("North American Breeding Bird Survey", dmeta2$Citation, ignore.case = T))
dmeta3=filter(dmeta2, !(ID %in% c(removed$ID, removed2$ID)))
#manually check remaining data!
write.csv(dmeta3, "data/LP/LP_datacheck.csv", row.names = F)
#removed additional model output and smoothed/interpolated survey data
dcheck=read.csv("data/LP/LP_datacheck_manual.csv")
keep=dcheck$ID[which(dcheck$exclude_manual=="no")]

dmeta3=filter(dmeta3, ID %in% keep)
dnestfilt=dnestfilt %>%
  filter(ID %in% dmeta3$ID)

#calculate monotonic trend
monotonic_eval=function(data) {
  cor(data$Year, data$Abundance, use="p", method="spearman")
}
dnestfilt$monotonic=map_dbl(dnestfilt$data_rescale, monotonic_eval)

# hightrend=filter(dnestfilt, abs(monotonic)>0.9)
# par(mfrow=c(4,3),mar=c(4,4,2,1))
# for(i in hightrend$ID) {
#   plot(Abundance~Year, data=subset(allseries, ID==i),
#        main=paste(i, dmeta$Common_name[dmeta$ID==i]), type="o", col="red")
# }
# hightrend2=filter(dnestfilt, abs(monotonic)<0.9 & abs(monotonic)>0.8)
# par(mfrow=c(4,3),mar=c(4,4,2,1))
# for(i in hightrend2$ID) {
#   plot(Abundance~Year, data=subset(allseries, ID==i),
#        main=paste(i, dmeta$Common_name[dmeta$ID==i]), type="o", col="orange")
# }
# decided not to remove these, will remove posthoc
# dnestfilt=filter(dnestfilt, abs(monotonic)<0.9)
# dmeta4=filter(dmeta3, ID %in% dnestfilt$ID)

#unnest tables
allseries_final=unnest(dnestfilt, data_rescale) %>% select(ID,Year:AbundanceLogScale)
allseries_meta=inner_join(dmeta, select(dnestfilt,ID,ndatapoints,data_rescale_case:monotonic)) %>%
  mutate(TaxonomicGroup=ifelse(Class=="Aves","Birds",ifelse(Class=="Mammalia","Mammals","Fishes")))

#add body size information to metadata
# #amniote traits (birds and mammals) from Myhrvold et al. 2025
# amtraits=read.csv("data/LP/lifehistory/Amniote_Database_Aug_2015.csv", na.strings = -999) %>%
#   mutate(Binomial2=paste(genus, species))
# #fish traits from Beukhof et al. 2019
# fish_lh=readxl::read_xlsx('data/LP/lifehistory/TraitCollectionFishNAtlanticNEPacificContShelf.xlsx')
#
# lp_lh=allseries_meta %>%
#   select(Binomial, Common_name, Family, Class) %>%
#   unique.data.frame() %>% arrange(Class, Binomial) %>%
#   mutate(Binomial2=sub("_"," ",Binomial))
# lp_lh=left_join(lp_lh, amtraits %>% select(Binomial2, Mass_g=adult_body_mass_g, MinAge_mo=female_maturity_d)) %>%
#   mutate(MinAge_mo=MinAge_mo/365*12)
#
# fish_lh=fish_lh %>% group_by(taxon) %>%
#   summarise(age.maturity=mean(age.maturity, na.rm=T),
#             length.max=mean(length.max, na.rm=T))
# lp_lh=lp_lh %>% left_join(fish_lh, by=c("Binomial2"="taxon"))
# lp_lh=lp_lh %>% mutate(MinAge_mo=ifelse(is.na(MinAge_mo), age.maturity*12, MinAge_mo))
# write.csv(lp_lh, "data/LP/LP_lifehistory.csv", row.names = F)

#manually edit (some species are in amniote database but names have changed), and add additional fish data
lp_lh=read.csv("data/LP/LP_lifehistory_manual.csv", stringsAsFactors = F)
allseries_meta=allseries_meta %>% left_join(lp_lh)
allseries_meta$SamplingInterval="annual"
allseries_meta$mo_timestep=12

#adjust columns for consistency
allseries_meta=allseries_meta %>%
  mutate(Id=paste('LP',ID, sep = "_"),Loc=substr(allseries_meta$Location,1,40), Name=Common_name) %>%
  relocate(Id,Loc,Name) %>%
  select(Id,Loc,Name,ID,Binomial,Common_name,Citation,SamplingInterval,mo_timestep,ndatapoints,longestrun,propmissing,propzeros,monotonic,
         Mass_g,MinAge_mo,TaxonomicGroup,Location,Country,Longitude,Latitude,data_rescale_case)
allseries_final=allseries_final %>% ungroup() %>%
  mutate(ID=paste('LP',ID, sep = "_")) %>%
  rename(Id=ID)

#write output
write.csv(allseries_final,"data/LP/LP_series.csv", row.names = F)
write.csv(allseries_meta,"data/LP/LP_seriesmetadata.csv", row.names = F)

table(allseries_meta$Class)
table(allseries_meta$TaxonomicGroup)

#plots ####
#allseries_meta=arrange(allseries_meta, TaxonomicGroup, Binomial)

# IDs=allseries_meta$ID
# par(mfrow=c(4,3),mar=c(4,4,2,1))
# for(i in IDs) {
#   plot(Abundance~Year, data=subset(allseries_final, ID==i),
#        main=paste(i, allseries_meta$Common_name[allseries_meta$ID==i]), type="o")
#   if(which(allseries_meta$ID==i)%%12==0) {
#     readline(prompt="Press [enter] to continue")
#   }
# }

#check for duplicates in GPDD ####

# exclude=c(27571,1047,6053,6098, 22059, 22062,22063,22073,22076)
#
# gpdd_d=read.csv("data/GPDD/gpdd_ts_metadata.csv")
# gpdd_d_ts=read.csv("data/GPDD/gpdd_timeseries.csv")
# gpdd_d_out=read.csv("results/gpdd_results.csv", stringsAsFactors = F)
# gpdd_d_sub=gpdd_d %>% filter(MainID %in% gpdd_d_out$MainID)
#
# allseries_meta=inner_join(dmeta, select(dnestfilt,ID,data_rescale_case:effectivenE2)) %>%
#   mutate(TaxonomicGroup=ifelse(Class=="Aves","Birds",ifelse(Class=="Mammalia","Mammals","Fishes")))
#
# allseries_meta$Binomial2=sub("_"," ",allseries_meta$Binomial)
# sharedsp=intersect(allseries_meta$Binomial2, gpdd_d_sub$TaxonName)
#
# test_gpdd=filter(gpdd_d_sub, TaxonName %in% sharedsp)
# test_lp=filter(allseries_meta, Binomial2 %in% sharedsp)
#
# #same species name
# #cross referenced location and citation for each
#
# #wolves in belarus
# #from the same paper, but values not identical, gpdd longer
# i=27571 #exclude
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=9393
# plot(Population~SampleYear, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
#
# #whooping cranes in aransas
# #different papers, same location, series have different year ranges
# i=1047 #exclude
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=1347
# plot(Population~SampleYear, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
#
# #sockeye salmon in skeena river
# #different papers, same location, LP series from diff tributaries, GPDD for whole river
# i=22059 #exclude
# i=22062
# i=22063
# i=22073
# i=22076
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=2759
# plot(Population~SampleYear, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
#
# #GPDD cited
#
# #sandgrouse in south africa
# #looks same, gpdd series is longer
# i=6053 #exclude
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=9793
# plot(Population~SampleYear, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
#
# #voles in finland (keep all)
# #LP has updated genus name. series do not look at all alike, diff site info
# #looking at original publication, it appears the series are indeed from 2 different sites
# i=5831
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=9920
# plot(Population~SeriesStep, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
#
# #pink salmon at bonneville dam
# #same, but GPDD series is longer
# i=6098 #exclude
# plot(Abundance~Year, data=subset(allseries_final, ID==i),
#      main=paste("LP", i, allseries_meta$Binomial2[allseries_meta$ID==i]), type="o")
# i=7067 #not in dataset, was omitted because reliability was -1
# plot(Population~SampleYear, data=subset(gpdd_d_ts, MainID==i),
#      main=paste("GPDD", i, gpdd_d_sub$TaxonName[gpdd_d_sub$MainID==i]), type="o")
