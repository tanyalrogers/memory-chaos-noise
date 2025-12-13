#plot results from all analyses

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
source("scripts/ggplot_themes.R")

#load data####

allresults_filt=read.csv("results/combined_results.csv")
#order levels
allresults_filt$JLEsign2=factor(allresults_filt$JLEsign2, levels = c("stable", "neutral", "unstable"))
allresults_filt$JLEsign2_mo=factor(allresults_filt$JLEsign2_mo, levels = c("stable", "neutral", "unstable"))

# marginal statistics ####

table(allresults_filt$Group)
table(allresults_filt$TaxonomicGroup)
table(allresults_filt$TaxonomicClass3)
table(allresults_filt$SamplingInterval)

#prop of series with E>1
allresults_filt %>% filter(E>1) %>% nrow() / nrow(allresults_filt)
#prop of series with theta>0
allresults_filt %>% filter(theta>0) %>% nrow() / nrow(allresults_filt)

#prop chaotic
table(allresults_filt$Group, allresults_filt$JLEsign)

#colors
stabilitycols=rev(RColorBrewer::brewer.pal(3, "Set1"))

# overall best R2 ####

#with smoother and quantreg, monthly
(f1a=ggplot(allresults_filt %>% filter(JLEmean_mo<1), aes(x=JLEmean_mo, y=smapoptR2)) +
    #facet_wrap(.~TaxonomicGroup, scales = "free_x") +
    geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
    #geom_smooth(color="black", method = "loess", span=1) +
    geom_quantile(color="black", formula = y ~ poly(x, 2), quantiles=c(0.8), lty=2) +
    geom_smooth(color="black", method="lm", formula = y ~ poly(x, 2), se = T, linewidth=0.75) +
    geom_point(aes(color=TaxonomicGroup, shape=TaxonomicGroup), alpha=0.6, size=2) +
    geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
    annotate("text", x=c(-0.4,0,0.4), y=1.05, label=c("stable","neutral","unstable"), color=stabilitycols) +
    scale_color_brewer(palette = "Dark2") +
    labs(x=expression(LE~(month^-1)), y="NSE (full model)", color="Taxonomic\nGroup", shape="Taxonomic\nGroup") +
    classic + removefacetbackground)
ggsave("figures/fig1.png", f1a, width = 6, height = 4)

# memory vs LE ####

(f1b=ggplot(allresults_filt %>% filter(JLEmean_mo<1), aes(x=JLEmean_mo, y=meanlag_st)) +
    geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
    geom_point(aes(fill=JLEsign2_mo), alpha=0.5, pch=21, size=2) +
    geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    labs(x=expression(LE~(month^-1)), y="Memory length (months)", fill="LE classification", color="E") +
    classic)
#ggsave("figures/meanlag_v_LE.png", f1b, width = 6, height = 4)

(f1c=ggplot(allresults_filt %>% filter(JLEmean_mo>0), aes(x=log(abs(1/JLEmean_mo)), y=log(meanlag_st))) +
    #facet_grid(.~JLEsign2_mo) +
    #geom_abline(slope = 1, intercept = 0, color="gray70") +
    geom_point(aes(color=TaxonomicGroup, shape=TaxonomicGroup),alpha=0.75, size=2) +
    geom_smooth(method = "lm", se=F, color="black", linewidth=0.75) +
    geom_hline(yintercept = 0, color="gray20") + # geom_vline(xintercept = 0, color="gray20") +
    scale_color_brewer(palette = "Dark2") +
    labs(x=expression(Lyapunov~horizon:~ln~1/LE~(month^-1)), y="ln Memory length (months)", color="Taxonomic\nGroup", shape="Taxonomic\nGroup") +
    classic)
#ggsave("figures/meanlag_v_LE_horizon.png", f1c, width = 6, height = 4)

f1b/f1c + plot_annotation(tag_levels = "a")
ggsave("figures/fig3.png", width = 6, height = 7)

# prediction improvement scatterplots ####

# * diff with best 1d model #####

# monthly, all datapoints, exclude outlier
p1s2=ggplot() +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(data=allresults_filt %>% filter(JLEmean_mo<1),
             aes(x=JLEmean_mo, y=diffbest1d, color=TaxonomicGroup, size=smapoptR2),
             alpha=0.4) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  annotate("text", x=c(-0.4,0,0.4), y=0.85, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(color="Taxonomic\nGroup", #fill="Sign.\nMemory",
       x=expression(LE~(month^-1)), y=expression(Delta~NSE~(full-"1d")), size="NSE (full)") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))
p1s2

#zoomed in
p1s3=ggplot(allresults_filt %>% filter(abs(JLEmean_mo)<0.2),
            aes(x=JLEmean_mo, y=diffbest1d, color=TaxonomicGroup, size=smapoptR2)) +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(alpha=0.4) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  annotate("text", x=c(-0.12,0,0.12), y=0.85, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(color="Taxonomic\nGroup", #fill="Sign.\nMemory",
       x=expression(LE~(month^-1)), y=expression(Delta~NSE~(full-"1d")), size="NSE (full)") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))
p1s3

# #colored by taxon
# p1=ggplot(allresults_filt %>% filter(E>1),
#        aes(x=JLEmean, y=diffbest1d, color=TaxonomicGroup, size=smapoptR2)) +
#   #facet_wrap(.~TaxonomicGroup, scales="free_x") +
#   #geom_smooth(color="black", method = "loess", span=0.9, show.legend = F) +
#   geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
#   geom_point(alpha=0.5) +
#   geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
#   scale_color_brewer(palette = "Set1") +
#   scale_size_continuous(range = c(1,5)) +
#   labs(color="Taxonomic\nGroup", #fill="Sign.\nMemory",
#        x=expression(LE~(timestep^-1)), y=expression(R[full]^2-R[1-d]^2), size=expression(R[full]^2)) +
#   classic + removefacetbackground +
#   guides(color = guide_legend(override.aes = list(size = 3)))
# p1
# ggsave("figures/diff1d_colortaxon.png", p1, width = 6, height = 4)
#
# #with monthly LE
# p1s=ggplot(allresults_filt %>% filter(E>1),
#           aes(x=JLEmean_mo, y=diffbest1d, color=SamplingInterval, size=smapoptR2)) +
#   #facet_wrap(.~TaxonomicGroup, scales="free_x") +
#   #geom_smooth(color="black", method = "loess", span=0.9, show.legend = F) +
#   geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
#   geom_point(alpha=0.5) +
#   geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
#   scale_color_brewer(palette = "Dark2") +
#   scale_size_continuous(range = c(1,5)) +
#   labs(color="Sampling\nInterval", #fill="Sign.\nMemory",
#        x=expression(LE~(month^-1)), y=expression(R[full]^2-R[1-d]^2), size=expression(R[full]^2)) +
#   classic + removefacetbackground +
#   guides(color = guide_legend(override.aes = list(size = 3)))
# p1s

#shape and color taxon, no size scaling
# ggplot(allresults_filt %>% filter(E>1),
#        aes(x=JLEmean, y=diffbest1d, color=TaxonomicGroup, shape=TaxonomicGroup)) +
#   #facet_wrap(.~TaxonomicGroup, scales="free_x") +
#   geom_point(alpha=0.75) +
#   geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
#   scale_color_brewer(palette = "Set1") +
#   labs(fill="Taxonomic\nGroup", #fill="Sign.\nMemory",
#        x=expression(LE~(timestep^-1)), y=expression(R[full]^2-R[1-d]^2), size=expression(R[full]^2)) +
#   classic + removefacetbackground +
#   guides(fill = guide_legend(override.aes = list(size = 3)))

# * diff with best AR(E) #####

# monthly, all datapoints, exclude outlier
p2s2=ggplot(allresults_filt %>% filter(JLEmean_mo<1),
            aes(x=JLEmean_mo, y=diffAR, color=TaxonomicGroup, size=smapoptR2)) +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(alpha=0.4) +
  #geom_smooth(color="black", method = "loess", span=1, show.legend = F) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  annotate("text", x=c(-0.4,0,0.4), y=0.85, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(color="Taxonomic\nGroup",
       x=expression(LE~(month^-1)), y=expression(Delta~NSE~(full-linear)), size="NSE (full)") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))
p2s2

#zoomed in
p2s3=ggplot(allresults_filt %>% filter(abs(JLEmean_mo)<0.2),
            aes(x=JLEmean_mo, y=diffAR, color=TaxonomicGroup, size=smapoptR2)) +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
  geom_point(alpha=0.4) +
  #geom_smooth(color="black", method = "loess", span=1, show.legend = F) +
  geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(range = c(1,5)) +
  annotate("text", x=c(-0.12,0,0.12), y=0.85, label=c("stable","neutral","unstable"), color=stabilitycols) +
  labs(color="Taxonomic\nGroup",
       x=expression(LE~(month^-1)), y=expression(Delta~NSE~(full-linear)), size="NSE (full)") +
  classic + removefacetbackground +
  guides(color = guide_legend(override.aes = list(size = 3)))
p2s3

# colored by taxon
# p2=ggplot(allresults_filt %>% filter(theta>0),
#        aes(x=JLEmean, y=diffAR, color=TaxonomicGroup, size=smapoptR2)) +
#   #facet_wrap(.~TaxonomicGroup, scales="free_x") +
#   geom_smooth(color="black", method = "loess", span=0.9, show.legend = F) +
#   geom_vline(xintercept = c(-0.02,0.02), color="gray70") +
#   geom_point(alpha=0.5) +
#   geom_hline(yintercept = 0, color="gray20") + geom_vline(xintercept = 0, color="gray20") +
#   scale_color_brewer(palette = "Set1") +
#   scale_size_continuous(range = c(1,5)) +
#   labs(color="Taxonomic\nGroup",
#        x=expression(LE~(timestep^-1)), y=expression(R[full]^2-R[AR]^2), size=expression(R[full]^2)) +
#   classic + removefacetbackground +
#   guides(color = guide_legend(override.aes = list(size = 3)))
# #ggsave("figures/diffAR_colortaxon.png", p3, width = 6, height = 4)

#boxplots ####

b10=ggplot(allresults_filt %>% filter(E>1), aes(fill=JLEsign2, y=diffbest1d, x=Group)) +
  #facet_wrap(.~Group, scales = "free_y") +
  geom_hline(yintercept = 0, color="gray20") +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_jitter(aes(group=JLEsign2), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0), alpha=0.3, color="grey40") +
  labs(fill="LE classification", x="Taxonomic Group", y=expression(Delta~NSE~(full-"1d"))) +
  classic + removefacetbackground

(b1=ggplot(allresults_filt, aes(fill=JLEsign2_mo, y=diffbest1d, x=Group)) +
  #facet_wrap(.~Group, scales = "free_y") +
  geom_hline(yintercept = 0, color="gray20") +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_jitter(aes(group=JLEsign2_mo), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0), alpha=0.3, color="grey40") +
  labs(fill="LE classification", x="Taxonomic Group", y=expression(Delta~NSE~(full-"1d"))) +
  classic + removefacetbackground)

b20=ggplot(allresults_filt %>% filter(theta>0), aes(fill=JLEsign2,y=diffAR, x=Group)) +
  #facet_wrap(.~Group, scales = "free_y") +
  geom_hline(yintercept = 0, color="gray20") +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_jitter(aes(group=JLEsign2), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), alpha=0.3, color="grey40") +
  labs(fill="LE classification", x="Taxonomic Group", y=expression(Delta~NSE~(full-linear))) +
  classic + removefacetbackground

(b2=ggplot(allresults_filt, aes(fill=JLEsign2_mo,y=diffAR, x=Group)) +
  #facet_wrap(.~Group, scales = "free_y") +
  geom_hline(yintercept = 0, color="gray20") +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_jitter(aes(group=JLEsign2_mo), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), alpha=0.3, color="grey40") +
  labs(fill="LE classification", x="Taxonomic Group", y=expression(Delta~NSE~(full-linear))) +
  classic + removefacetbackground)

#combo plots ####

f31=(p1s3+p2s3) + plot_layout(guides = "collect")
f32=(b1+b2) + plot_layout(guides = "collect")
f31/f32 + plot_annotation(tag_levels = list(c("a","c","b","d")))
ggsave("figures/fig2.png", width = 9, height = 7)

(p1s2+p2s2) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")
ggsave("figures/dNSE_fullrange.png", width = 9, height = 4)

b10+b20 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
ggsave("figures/dNSE_boxplot_omiting.png", width = 9, height = 4)

#LE and body size ####
(bs1=ggplot(allresults_filt %>% filter(JLEmean_mo>0 & !is.na(Mass_g)),
       aes(y=log10(JLEmean_mo), x=log10(Mass_g))) +
  ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
  geom_smooth(method="lm", se = F, color="black", linewidth=0.75) +
  geom_point(aes(color=TaxonomicGroup, shape=TaxonomicGroup), alpha=0.75, size=2) +
  classic + labs(color="Taxonomic\nGroup", shape="Taxonomic\nGroup") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=10)))

#by source
# ggplot(allresults_filt %>% filter(JLEmean_mo>0 & !is.na(Mass_g)),
#        aes(y=log10(JLEmean_mo), x=log10(Mass_g))) +
#   ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
#   geom_smooth(method="lm", se = F, color="black") +
#   geom_point(aes(fill=Source), color="black", alpha=0.75, pch=21, size=3.5) +
#   classic + labs(fill="Source") + scale_fill_manual(values = c("firebrick","cornflowerblue","gold","forestgreen")) +
#   # guides(fill = guide_legend(override.aes = list(size = 3.5, shape=21))) +
#   # guides(shape = guide_legend(override.aes = list(stroke = c(0.5,0.5,1.25)))) +
#   theme(axis.title = element_text(size=12), axis.text = element_text(size=10))

#Memory and body size ####
(bs2=ggplot(allresults_filt %>% filter(!is.na(Mass_g) & meanlag_st>12),
            aes(y=meanlag_st, x=log10(Mass_g), group=JLEsign2_mo)) +
   #facet_grid(.~JLEsign2_mo)+
   ylab("Memory length (months)") +
   xlab(expression(~log[10]~Mass~(g))) +
   geom_smooth(aes(color=JLEsign2_mo), method="lm", se = F, linewidth=0.75) +
   geom_point(aes(color=JLEsign2_mo, group=JLEsign2_mo), alpha=0.75, size=2) +
   classic + labs(color="LE classification", shape="Taxonomic\nGroup") +
   scale_color_brewer(palette = "Set1", direction = -1) +
   theme(axis.title = element_text(size=12), axis.text = element_text(size=10)))

bs1/bs2 + plot_annotation(tag_levels = "a") #+ plot_layout(guides = "collect")
ggsave("figures/fig4.png", width = 6, height = 6)

mm=lm(meanlag_st~log10(Mass_g), data=allresults_filt %>% filter(!is.na(Mass_g) & meanlag_st>12))
anova(mm)

allresults_filt %>% filter(!is.na(Mass_g)) %>% nrow()
allresults_filt %>% nrow()

#Supplemental Figs ####

# histograms of E, theta ####
(h1=ggplot(allresults_filt,
          aes(x=factor(theta), fill = JLEsign2_mo)) +
  #facet_grid(JLEsign2_mo~SamplingInterval2, scales="free_y") +
  facet_grid(JLEsign2_mo~., scales="free_y") +
  #geom_vline(xintercept = 0, color="gray70") +
  geom_bar(show.legend = F) + scale_fill_manual(values = stabilitycols) +
  labs(fill=expression(Smap~theta), #fill="Sign.\nMemory",
       x=expression(Nonlinearity~(theta))) +
  classic + removefacetbackground + scale_y_continuous(expand = expansion(mult=c(0,0.1))))

(h2=ggplot(allresults_filt,
          aes(x=E, fill = JLEsign2_mo)) +
  #facet_grid(JLEsign2_mo~SamplingInterval2, scales="free_y") +
  facet_grid(JLEsign2_mo~., scales="free_y") +
  geom_bar(show.legend = F) + scale_fill_manual(values = stabilitycols) +
  labs(x=expression(Dimensionality~(E))) +
  classic + removefacetbackground +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  scale_x_continuous(breaks = 1:6))

h2+h1 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
ggsave("figures/histo_combined.png", width = 9, height = 5)

with(allresults_filt, table(JLEsign2_mo, E))
with(allresults_filt, table(JLEsign2_mo, theta))
with(allresults_filt, table(JLEsign2_mo))

# R2 smap vs. R2 for 1-d and AR models ####

mc1=ggplot(allresults_filt, aes(y=smapoptR2, x=smap1dR2, color=TaxonomicGroup, shape=TaxonomicGroup)) +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_point(alpha=0.75, size=3) +
  geom_abline(intercept = 0, slope=1, color="gray20", lwd=1) +
  scale_color_brewer(palette = "Dark2") +
  labs(color="Taxonomic\nGroup",shape="Taxonomic\nGroup",
       y="NSE (full)", x="NSE (1d)") +
  classic + removefacetbackground +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

mc2=ggplot(allresults_filt, aes(y=smapoptR2, x=ARmodR2, color=TaxonomicGroup, shape=TaxonomicGroup)) +
  #facet_wrap(.~TaxonomicGroup, scales="free_x") +
  geom_point(alpha=0.75, size=3) +
  geom_abline(intercept = 0, slope=1, color="gray20", lwd=1) +
  scale_color_brewer(palette = "Dark2") +
  labs(color="Taxonomic\nGroup",shape="Taxonomic\nGroup",
       y="NSE (full)", x="NSE (linear)") +
  classic + removefacetbackground +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

mc1+mc2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
ggsave("figures/NSE_compare.png", width = 9, height = 4)
