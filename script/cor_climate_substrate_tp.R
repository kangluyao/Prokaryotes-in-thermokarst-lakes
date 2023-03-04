# test the correlation between the climate and substrate with TP dataset
metadata %>% filter(Region == 'Tibetan Plateau') %>%
  dplyr::select(BIX) %>%
  dplyr::summarise_all(list(min = min, max = max, 
                            media = median, mean = mean, 
                            sd = sd, se = stderr), na.rm = TRUE)

cor_map_bix <- ggplot(metadata, aes(x = MAP, y = BIX)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  ggpubr::stat_cor(aes(MAP, BIX, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) + 
  xlab('MAP (mm)')+ylab('BIX') +
  scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

cor_map_bix_den <- ggExtra::ggMarginal(cor_map_bix, type = "histogram", 
           margins = "y", fill= 'tomato3')

cor_map_hix <- ggplot(metadata, aes(x = MAP, y = HIX)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  ggpubr::stat_cor(aes(MAP, HIX, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", 
                   label.x.npc = 0.1, label.y.npc = 0.1, size = 5) + 
  xlab('MAP (mm)')+ylab('HIX') +
  scale_y_continuous(expand = c(0,0), limits = c(0.3, 1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

cor_map_hix_den <- ggExtra::ggMarginal(cor_map_hix, type = "histogram", 
                                       margins = "y", fill= 'tomato3')

library(cowplot)
cor_climate_dom <- ggdraw() +
  draw_plot(cor_map_bix_den, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(cor_map_hix_den, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0, 0.5), y = c(1, 1))
cor_climate_dom
