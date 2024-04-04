plot.pcaxis <- function(axis = 1)
{
  
  pc.ax <- pc.coord[,paste("Axis",axis,sep="")]
  ax <- as.vector(pc.ax[[1]])
  names(ax) <-pc.coord$label
  fit.ax <- phytools::fastAnc(acou.tree, ax, vars=TRUE, CI=TRUE)
  
  td.ax <- data.frame(node = nodeid(acou.tree, names(ax)),
                      trait = ax)
  nd.ax <- data.frame(node = names(fit.ax$ace), trait = fit.ax$ace)
  
  d.ax <- rbind(td.ax, nd.ax)
  d.ax$node <- as.numeric(d.ax$node)
  tree.ax <- full_join(acou.tree, d.ax, by = 'node')
  
  p1 <- ggtree(tree.ax, aes(color=trait), layout = 'circular', 
               ladderize = FALSE, continuous = 'colour', size=1)+
    scale_color_viridis_c()+
    geom_tiplab(size=2, aes(angle=angle))+
    labs(title = paste("Acoustic trait space - axis",axis,sep=" "))
  
  print(p1)
  
}
