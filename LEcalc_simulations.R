# Applies Jacobian LE chaos detection method to simulated series

library(rEDM)
library(dplyr)
library(tidyr)
library(purrr)

source("scripts/chaos_localstability_methods.R")

ts_d=read.csv("./data/simulations/NoisyRM.csv", header = F)
ts_d=exp(ts_d)
ts_d=ts_d %>% mutate(Time=1:n()) %>% pivot_longer(names_to = "ID", values_to = "Abundance", cols = 1:60) %>%  group_by(ID) %>% nest(.key="data_rescale") %>%
  mutate(data_rescale=map(data_rescale, as.data.frame))
#for updating
#dnestfilt=ts_d

ts_results=select(ts_d, ID)

#### Jacobian LE Method ####

#fit models
#data are already log transformed
ts_results$modelresultsbest=map(ts_d$data_rescale, smap_model_options, y="Abundance", model=5, seaspred=F, taufix=1) #gr-log
ts_d$modelform=map_chr(ts_results$modelresultsbest, ~.x$form)

#pull out R2 values
ts_d$bestR2m=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$R2model)
ts_d$bestR2=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$R2abund)

#get E, tau, theta values from model
ts_d$E=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$E)
ts_d$tau=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$tau)
ts_d$theta=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$theta)

#get jacobian matrices
ts_results$jacobians=map(ts_results$modelresultsbest, getJacobians)

#get global LE
ts_results$GlobalLE=map2(ts_results$modelresultsbest, ts_results$jacobians, LEshift)
ts_d$JLEmin=map_dbl(ts_results$GlobalLE, ~.x$minci) #LE lower confidence bound (*this is the LE estimate to use!*)
ts_d$JLEmean=map_dbl(ts_results$GlobalLE, ~.x$minmean) #LE mean (*this is the LE estimate to use!*)
ts_d$JLEsign=ifelse(ts_d$JLEmin>0.01, "chaotic", "not chaotic")

#### Fit alternative models ####

#1d smap
ts_results$smap1d=map(ts_d$data_rescale, smap_model_options, Efix=1, model=5, taufix=1, y="Abundance")
ts_d$smap1dR2=map_dbl(ts_results$smap1d, ~.x$modelstats$R2abund)

#best AR
ts_results$ARmod=map(ts_d$data_rescale, smap_model_options, thetafix=0, model=5, taufix=1, y="Abundance")
ts_d$ARmodR2=map_dbl(ts_results$ARmod, ~.x$modelstats$R2abund)

#### Export Results ####

#save results
save(ts_d, ts_results, file = "results/simulations_LE_12-16-25.Rdata")

#global LE and other site-level values
exportres1=ts_d %>% select(ID, R2abund=bestR2, R2gr=bestR2m, modelform,
                           E, tau, theta, JLEmin, JLEmean, JLEsign,
                           smap1dR2, smap1dR2, ARmodR2)
write.csv(exportres1, "results/simulations_LEresults.csv", row.names = F)

#### Plots ####

library(ggplot2)
library(patchwork)
source("scripts/ggplot_themes.R")
stabilitycols=rev(RColorBrewer::brewer.pal(3, "Set1"))

exportres1$diffbest1d=exportres1$R2abund-exportres1$smap1dR2
exportres1$diffAR=exportres1$R2abund-exportres1$ARmodR2
exportres1$obserror=rep(c(0,0.05,0.1), each=20)

#diff NSE 1d
p1=ggplot() +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(data=exportres1 %>% filter(R2abund>0.15),
             aes(x=JLEmean, y=diffbest1d, size=R2abund, color=factor(obserror)),
             alpha=0.4) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  annotate("text", x=c(-0.2,0,0.15), y=0.63, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(x=expression(LE~(timestep^-1)), y=expression(Delta~NSE~(full-"1d")), size="NSE (full)", color="Obs. error") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))

#diff NSE AR
p2=ggplot() +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(data=exportres1 %>% filter(R2abund>0.15),
             aes(x=JLEmean, y=diffAR, size=R2abund, color=factor(obserror)),
             alpha=0.4) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  annotate("text", x=c(-0.2,0,0.15), y=1.13, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(x=expression(LE~(timestep^-1)), y=expression(Delta~NSE~(full-linear)), size="NSE (full)", color="Obs. error") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))

#full model NSE
p3=ggplot() +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(data=exportres1 %>% filter(R2abund>0.15),
             aes(x=JLEmean, y=R2abund, color=factor(obserror)),
             alpha=0.4, size=3) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  annotate("text", x=c(-0.2,0,0.15), y=1.03, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(x=expression(LE~(timestep^-1)), y="NSE (full)", color="Obs. error") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))

#linear vs full NSE
p4=ggplot() +
  geom_abline(intercept = 0, slope=1, color="gray20", lwd=1) +
  geom_point(data=exportres1 %>% filter(R2abund>0.15),
             aes(x=ARmodR2, y=R2abund, color=factor(obserror)),
             alpha=0.4, size=3) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  labs(x="NSE (linear)", y="NSE (full)", color="Obs. error") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))

(p4+p3)/(p1+p2)+ plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
ggsave("figures/sim_results.png", width = 8, height = 6.5)
