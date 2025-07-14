#InsectChange time series - data processing

library(dplyr)
library(tidyr)

#downloaded data files
dsources=read.csv("C:/Users/trogers/Documents/Nonstationarity/InsectChange/DataSources.csv")
dplot=read.csv("C:/Users/trogers/Documents/Nonstationarity/InsectChange/PlotData.csv")
dtimeseries=read.csv("C:/Users/trogers/Documents/Nonstationarity/InsectChange/InsectAbundanceBiomassData.csv")

un4=function(m,x) { #number of months sampled at least 4 times
  if(length(unique(m>1))) {
    ts=!is.na(x)
    tab=table(m,ts)
    return(length(which(tab[,"TRUE"]>4)))
  } else {
    return(1)
  }
}

longestrun=function(x) { #longest string of nonmissing values
  len=rle(!is.na(x))
  len2=len$lengths[len$values==TRUE]
  return(max(len2))
}

dlength=dtimeseries %>% group_by(DataSource_ID, Plot_ID) %>%
  summarise(nvalues=n(),ndatapoints=length(which(!is.na(Number))),
            uniqueperiods=length(unique(Period)),
            uniqueperiods4=un4(Period,Number),
            uniqueyears=length(unique(Year)), .groups = "drop") %>%
  mutate(SamplingIntervalOriginal=ifelse(uniqueperiods==1,"annual","monthly"),
         propmissing=(nvalues-ndatapoints)/nvalues)
#less than 22% missing, either 40 years or least 10 months sampled at least 4 times
dlengthlong=subset(dlength, (uniqueyears>=40 | uniqueperiods4>=10) & propmissing<0.22)
dsourceslong=subset(dsources, DataSource_ID %in% dlengthlong$DataSource_ID) %>% left_join(dlengthlong) %>%
  subset(Plot_ID %in% dlengthlong$Plot_ID)
dplotlong=subset(dplot, Plot_ID %in% dlengthlong$Plot_ID) %>%
  left_join(dsourceslong)

insectseries=subset(dtimeseries, Plot_ID %in% dlengthlong$Plot_ID)
unique(insectseries$Period)

#annual series
annualseries=subset(insectseries, Plot_ID %in% dlengthlong$Plot_ID[dlengthlong$uniqueperiods==1])
annualseries=annualseries %>% mutate(Id=paste(DataSource_ID, Plot_ID, Stratum, sep="_"), Month=1)
#there are no missing years
#check for 40 contiguous
annlongestrun=annualseries %>% group_by(Id) %>% summarize(longestrun=longestrun(Number))
annkeep=subset(annlongestrun, longestrun>=40)
annualseries2=subset(annualseries, Id %in% annkeep$Id)
#final series
annualseries_final=select(annualseries2, Id, DataSource_ID, Plot_ID, Stratum, Year, Month, Abundance=Number)
annualseries_un=unique.data.frame(select(annualseries_final, Id, DataSource_ID, Plot_ID, Stratum)) %>% na.omit() %>%
  mutate(SamplingInterval="annual") %>% left_join(annkeep)

#monthly series
monthlyseries=subset(insectseries, Plot_ID %in% dlengthlong$Plot_ID[dlengthlong$uniqueperiods>=10])
monthlyseries=monthlyseries %>% mutate(Id=paste(DataSource_ID, Plot_ID, Stratum, sep="_"))
monthlyseries$Period=gsub(" ", "", monthlyseries$Period)
monthlyseries$Month=as.numeric(ifelse(monthlyseries$Period %in% month.name, match(monthlyseries$Period,month.name),
                                      ifelse(monthlyseries$Period %in% tolower(month.name), match(monthlyseries$Period,tolower(month.name)),
                                             monthlyseries$Period)))
monthlyseries=monthlyseries %>% arrange(DataSource_ID, Plot_ID, Stratum, Year, Month) %>%
  filter(MetricAB=="abundance") #exclude biomass (duplicates for each time step in some series)
#average multiple measurements per month
monthlyseries2=monthlyseries %>% group_by(Id, DataSource_ID, Plot_ID, Stratum, Year, Month) %>%
  summarise(Number=mean(Number, na.rm = T), .groups = "drop")
#fill in missing months
tll=split(monthlyseries2, monthlyseries2$Id)
tllfill=lapply(tll, FUN=function(d) {
  year.range=range(d$Year)
  dates=expand.grid(Year=year.range[1]:year.range[2], Month=1:12)
  dfill=full_join(dates, d, by = c("Year", "Month")) %>% arrange(Year, Month)
  dfill=dfill[(first(which(!is.na(dfill$Id)))):(last(which(!is.na(dfill$Id)))),]
  dfill$Id=as.character(na.omit(unique(dfill$Id)))
  return(dfill)
})
monthlyseriesfill=bind_rows(tllfill) %>% fill(DataSource_ID, Plot_ID, Stratum)
#remove series with too much missing data
monthlyseriessum=monthlyseriesfill %>%
  group_by(Id) %>%
  summarise(propmissing=length(which(is.na(Number)))/n(), .groups = "drop") %>%
  filter(propmissing<0.22)
monthlyseries2=monthlyseriesfill %>% filter(Id %in% monthlyseriessum$Id)
#check for 40 contigouous
monlongestrun=monthlyseries2 %>% group_by(Id) %>% summarize(longestrun=longestrun(Number))
monkeep=subset(monlongestrun, longestrun>=40)
monthlyseries3=subset(monthlyseries2, Id %in% monkeep$Id)
#final series
monthlyseries_final= monthlyseries3 %>% select(Id, DataSource_ID, Plot_ID, Stratum, Year, Month, Abundance=Number)
monthlyseries_un=unique.data.frame(select(monthlyseries_final, Id, DataSource_ID, Plot_ID, Stratum)) %>% na.omit() %>%
  mutate(SamplingInterval="monthly") %>% left_join(monkeep)

#annualized series
annualizedmonthlyseries=subset(insectseries, Plot_ID %in% dlengthlong$Plot_ID[dlengthlong$uniqueperiods<10 & dlengthlong$uniqueperiods>1])
annualizedmonthlyseries=annualizedmonthlyseries %>% mutate(Id=paste(DataSource_ID, Plot_ID, Stratum, sep="_"))
annualizedmonthlyseries2=annualizedmonthlyseries %>% group_by(Id, DataSource_ID, Plot_ID, Stratum, Year) %>%
  summarise(Number=mean(Number, na.rm = T), .groups = "drop") %>%
  mutate(Month=1)
#there are no missing years
#check for 40 contigouous
annmonlongestrun=annualizedmonthlyseries2 %>% group_by(Id) %>% summarize(longestrun=longestrun(Number))
annmonkeep=subset(annmonlongestrun, longestrun>=40)
annualizedmonthlyseries3=subset(annualizedmonthlyseries2, Id %in% annmonkeep$Id)
#final series
annualizedmonthlyseries_final=select(annualizedmonthlyseries3, Id, DataSource_ID, Plot_ID, Stratum, Year, Month, Abundance=Number)
annualizedmonthlyseries_un=unique.data.frame(select(annualizedmonthlyseries_final, Id, DataSource_ID, Plot_ID, Stratum)) %>% na.omit() %>%
  mutate(SamplingInterval="annual") %>% left_join(annmonkeep)

#final series timeseries
allseries_final=rbind(annualseries_final, annualizedmonthlyseries_final, monthlyseries_final)
#rescale
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
rescale=function(x, log=T, diag=F) {
  mind=min(x, na.rm=T)
  minnz=min(x[x>0], na.rm=T)
  if(mind>0) {
    xScale=x/sd(x, na.rm=T)
    case=1 #no zeros, nothing added
  } else {
    if(all(is.wholenumber(x), na.rm=T) & minnz<100) {
      xScale=(x+1)/sd(x, na.rm=T)
      case=2 #zeros, integers, plus 1, unless minnz is >100
    } else {
      xScale=(x+minnz)/sd(x, na.rm=T)
      case=3 #zeros, non-integers, plus min non-zero value
    }
  }
  if(log) {
    xScale=log(xScale)
  }
  if(diag) return(case)
  else return(xScale)
}
allseries_final=allseries_final %>% group_by(Id) %>%
  mutate(AbundanceScale=rescale(Abundance, log=F),
         AbundanceLogScale=rescale(Abundance, log=T))
rescalecases=allseries_final %>% group_by(Id) %>%
  summarise(data_rescale_case=rescale(Abundance, diag=T))
monotonic=allseries_final %>% group_by(Id) %>%
  mutate(decYear=Year+(Month-1)/12) %>%
  summarise(monotonic=cor(decYear, Abundance, use="p", method="spearman"))
#final series metadata
allseries_meta=rbind(annualseries_un, annualizedmonthlyseries_un, monthlyseries_un) %>%
  mutate(mo_timestep=ifelse(SamplingInterval=="monthly",1,12)) %>%
  left_join(rescalecases) %>% left_join(monotonic) %>%
  left_join(dsourceslong) %>% left_join(dplotlong %>% select(1:8,14:16))

#add body mass estimates
allseries_meta$Mass_g=ifelse(allseries_meta$InvertebrateGroup=="Mosquitos",0.003,
                             ifelse(allseries_meta$InvertebrateGroup=="Butterflies",0.5,NA))
allseries_meta$MinAge_mo=NA

#adjust cols for consistency
allseries_final=select(allseries_final, Id, Year, Month, Abundance, AbundanceScale, AbundanceLogScale) %>% ungroup()
allseries_meta=allseries_meta %>% mutate(Name=InvertebrateGroup, Loc=Location) %>%
  relocate(Id,Loc,Name) %>%
  select(Id,Loc,Name,DataSource_ID,Plot_ID,Stratum,InvertebrateGroup,DataSourceName,Reference,OpenAccessLicense,SamplingInterval,mo_timestep,ndatapoints,propmissing,longestrun,monotonic,
         Mass_g,MinAge_mo,Location,NationState,Longitude,Latitude,data_rescale_case)

#export data
write.csv(allseries_final, "data/insects/insects_series.csv", row.names = F)
write.csv(allseries_meta, "data/insects/insects_seriesmetadata.csv", row.names = F)

# make some plots

# IDs=allseries_meta$Id
# allseries_final$decYear=allseries_final$Year+(allseries_final$Month-1)/12
# par(mfrow=c(3,2),mar=c(4,4,2,1))
# for(i in IDs) {
#   plot(log(AbundanceScale)~decYear, data=subset(allseries_final, Id==i), main=i, type="o")
# }
# par(mfrow=c(3,2),mar=c(4,4,2,1))
# for(i in IDs) {
#   plot(AbundanceScale~decYear, data=subset(allseries_final, Id==i), main=i, type="o")
# }
# unique(allseries_meta$DataSource_ID)
