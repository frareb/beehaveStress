rm(list=ls())
dev.off()

library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

load("./input_data/out.RData")

heatmap_pop_plot <- out[!is.na(out$AUCpop),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUCpop))+
  ggtitle('Population size')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluepop<.05, pvaluepop>.01), aes(label = round(AUCpop, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluepop<.01), aes(label = round(AUCpop, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(0.5, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_pop_plot
save(heatmap_pop_plot, file = "New_graphs/Figure_pop_v1.RData")


heatmap_larve_plot <- out[!is.na(out$AUClarve),] %>% 
  ggplot(aes(x = ttime_n, y = interv, fill = AUClarve))+
  ggtitle('Number of larvae')+
  geom_tile()+
  geom_label(data = . %>% filter(pvaluelarve<.05, pvaluelarve>.01), aes(label = round(AUClarve, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pvaluelarve<.01), aes(label = round(AUClarve, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = "RdBu", limits = c(0.5, 1), oob = scales::squish, name = 'AUC')+
  scale_x_continuous(breaks = c(120, 151, 181, 212), labels = c('Ap', "My", "Jn", "Jl"), name = "Day of observation", sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(breaks = seq(50, 250, length.out=5), name = 'Prediction (Day)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = 'top', legend.justification = c(1,1))
heatmap_larve_plot
save(heatmap_larve_plot, file = "New_graphs/Figure_larve_v1.RData")

heatmap_pop_larve_plot <- plot_grid(heatmap_pop_plot + theme(axis.text.y.right = element_blank(), axis.title.y.right = element_blank()),
                                    heatmap_larve_plot + theme(axis.text.y.left = element_blank(), axis.title.y.left = element_blank()), 
                                    ncol = 2)
heatmap_pop_larve_plot
save(heatmap_pop_larve_plot, file = "New_graphs/Figure_pop_larve_v1.RData")



load("./New_graphs/Figure_pop_v1.RData")
load("./New_graphs/Figure_larve_v1.RData")
load("./New_graphs/Figure_pop_larve_v1.RData")

pdf(file = "./New_graphs/heatmap.pdf", width = 12)
  heatmap_pop_plot
  heatmap_larve_plot
  heatmap_pop_larve_plot
dev.off()