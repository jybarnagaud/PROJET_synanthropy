library("ggplot2")
library("magrittr")
library("dplyr")

x = get(load("outputs/SSI_CARTNAT_LAYER_4"))

effsize_res <- x$effSizes
sp_ssi <- x$speciesScores

sub_effsize_res <- effsize_res %>% 
  filter(Resolution == 1) %>% 
  left_join(sp_ssi, by = c("Species", "Resolution"))

# make the score a factor
sub_effsize_res$Index <- as.factor(sub_effsize_res$Index)

# plot the results
ggplot(sub_effsize_res, 
       aes(x = reorder(Species, -effsize), 
           y = -effsize, fill = Index)) +
  geom_hline(yintercept = 0.0, color = "darkgrey", 
             linewidth = 0.8, linetype = "dashed") +
  geom_boxplot() + 
  coord_flip() +
  scale_fill_brewer(name = "SSI", palette = "RdYlGn") +
  ylab("Effect size") +
  xlab("Species") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = c(0.9, 0.2)) +
  theme_bw()
