rm(list=ls())
# dev.off()

library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

load("./input_data/out.RData")

#---------------------------------------------------------------------
# Heatmaps adapted from Cabon et al. 2022
#---------------------------------------------------------------------

limcol <- min( c(out$AUCpop, 
                 out$AUClarve, 
                 out$AUChoneyenergy, 
                 out$AUCAFF, 
                 out$AUCLPS, 
                 out$AUCtotevendtoday), 
               na.rm=T)

#---------------------------------------------------------------------
# Population size + Number of larvae

# Population size
heatmap_pop_plot <- out[!is.na(out$AUCpop),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCpop))+
  ggtitle('Population size')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluepop<.05, pvaluepop>.01), aes(label = round(AUCpop, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluepop<.01), aes(label = round(AUCpop, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_pop_plot
save(heatmap_pop_plot, file = "New_graphs/Figure_pop_v1.RData")

# Number of larvae
heatmap_larve_plot <- out[!is.na(out$AUClarve),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUClarve))+
  ggtitle('Number of larvae')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluelarve<.05, pvaluelarve>.01), aes(label = round(AUClarve, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluelarve<.01), aes(label = round(AUClarve, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = "RdBu", limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_larve_plot
save(heatmap_larve_plot, file = "New_graphs/Figure_larve_v1.RData")

# Population size + Number of larvae
heatmap_pop_larve_plot <- plot_grid(heatmap_pop_plot + theme(axis.text.y.right = element_blank(), axis.title.y.right = element_blank()),
                                    heatmap_larve_plot + theme(axis.text.y.left = element_blank(), axis.title.y.left = element_blank()), 
                                    ncol = 2)
heatmap_pop_larve_plot
save(heatmap_pop_larve_plot, file = "New_graphs/Figure_pop_larve_v1.RData")

#---------------------------------------------------------------------
# Additional indicators

# Total weight
heatmap_weight_plot <- out[!is.na(out$AUCweight),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCweight))+
  ggtitle('Total weight')+
  geom_tile()+
  geom_label(data = . %>% filter(pvalueweight<.05, pvalueweight>.01), aes(label = round(AUCweight, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvalueweight<.01), aes(label = round(AUCweight, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_weight_plot
save(heatmap_weight_plot, file = "New_graphs/Figure_weight_v1.RData")

# Honey store
heatmap_honey_plot <- out[!is.na(out$AUChoneyenergy),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUChoneyenergy))+
  ggtitle('Honey store')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluehoneyenergy<.05, pvaluehoneyenergy>.01), aes(label = round(AUChoneyenergy, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluehoneyenergy<.01), aes(label = round(AUChoneyenergy, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_honey_plot
save(heatmap_honey_plot, file = "New_graphs/Figure_honey_v1.RData")

# Age of first foraging
heatmap_AFF_plot <- out[!is.na(out$AUCAFF),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCAFF))+
  ggtitle('Age of first foraging')+
  geom_tile()+
  geom_label(data = . %>% filter(pvalueAFF<.05, pvalueAFF>.01), aes(label = round(AUCAFF, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvalueAFF<.01), aes(label = round(AUCAFF, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_AFF_plot
save(heatmap_AFF_plot, file = "New_graphs/Figure_AFF_v1.RData")

# Life span period
heatmap_LSP_plot <- out[!is.na(out$AUCLPS),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCLPS))+
  ggtitle('Life span period')+
  geom_tile()+
  geom_label(data = . %>% filter(pvalueLSP<.05, pvalueLSP>.01), aes(label = round(AUCLPS, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvalueLSP<.01), aes(label = round(AUCLPS, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_LSP_plot
save(heatmap_LSP_plot, file = "New_graphs/Figure_LSP_v1.RData")

# Daily flight activity
heatmap_flight_plot <- out[!is.na(out$AUCtotevendtoday),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCtotevendtoday))+
  ggtitle('Daily flight activity')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluetotevendtoday<.05, pvaluetotevendtoday>.01), aes(label = round(AUCtotevendtoday, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluetotevendtoday<.01), aes(label = round(AUCtotevendtoday, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_flight_plot
save(heatmap_flight_plot, file = "New_graphs/Figure_flight_v1.RData")

#---------------------------------------------------------------------
# Heatmaps with signe data
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Population size + Number of larvae

out$signepop_char <- factor(NA, levels = c('positive', 'negative'))
out$signelarve_char <- factor(NA, levels = c('positive', 'negative'))
out$signepop_char <- ifelse(out$signepop==1, 'positive', 'negative')
out$signelarve_char <- ifelse(out$signelarve==1, 'positive', 'negative')

# Population size
heatmap_pop_plot_v2 <- out[!is.na(out$AUCpop),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCpop))+
  ggtitle('Population size')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluepop<.05, pvaluepop>.01), aes(label = round(AUCpop, 2) %>% str_replace("0",""), colour=signepop_char), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluepop<.01), aes(label = round(AUCpop, 2) %>% str_replace("0",""), colour=signepop_char), label.size = 0, size=3, fill = "transparent", fontface='bold')+
  scale_colour_manual(values=c('positive'='deeppink', 'negative'='blue'), name='Effect\nsign' ) +
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_pop_plot_v2
save(heatmap_pop_plot_v2, file = "New_graphs/Figure_pop_v2.RData")

# Number of larvae
heatmap_larve_plot_v2 <- out[!is.na(out$AUClarve),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUClarve))+
  ggtitle('Number of larvae')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluelarve<.05, pvaluelarve>.01), aes(label = round(AUClarve, 2) %>% str_replace("0",""), colour=signelarve_char), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluelarve<.01), aes(label = round(AUClarve, 2) %>% str_replace("0",""), colour=signelarve_char), label.size = 0, size=3, fill = "transparent", fontface='bold')+
  scale_colour_manual(values=c('positive'='deeppink', 'negative'='blue'), name='Effect\nsign' ) +
  coord_cartesian(expand = F)+
  scale_fill_distiller(breaks = c(0.5, 0.75, 1.0), type ="div", palette = 4, limits = c(limcol, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_larve_plot_v2
save(heatmap_larve_plot_v2, file = "New_graphs/Figure_larve_v2.RData")


#---------------------------------------------------------------------
# PDF Exportation
#---------------------------------------------------------------------

load("./New_graphs/Figure_pop_v1.RData")
load("./New_graphs/Figure_larve_v1.RData")
load("./New_graphs/Figure_pop_larve_v1.RData")

pdf(file = "./New_graphs/heatmap.pdf", width = 6)
  heatmap_pop_plot
  heatmap_pop_plot_v2
  heatmap_larve_plot
  heatmap_larve_plot_v2
  # heatmap_pop_larve_plot
  heatmap_weight_plot
  heatmap_honey_plot
  heatmap_AFF_plot
  heatmap_LSP_plot
  heatmap_flight_plot
dev.off()
