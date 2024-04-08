library(mapview)
library(ade4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(PerformanceAnalytics)

stoc <- read.csv("data/stoc_for_rql.csv")
stoc.eps <-
  stoc[, c("birdlife_sci_name", "id_carre_num",  "X", "Y")]

#--------------------------------#
#### check nb plots / species ####
#--------------------------------#

car.sp <- unique(stoc.eps[,c("birdlife_sci_name","id_carre_num")])
n.car.sp <- tapply(car.sp$id_carre,car.sp$birdlife_sci_name,"length")
sort(n.car.sp)

tot.car <- length(unique(stoc.eps$id_carre_num))
prop.sp.car <- 100*n.car.sp / tot.car
sort(prop.sp.car)

#-------------#
#### plots ####
#-------------#

stoc.id <- unique(stoc.eps[, c("id_carre_num", "X", "Y")])

mapview(stoc.id,
        xcol = 'X',
        ycol = 'Y',
        crs = 'epsg:27572')

#-----------------#
#### landscape ####
#-----------------#

land.keep <- c(
  "id_carre_num",
  "cartnat_layer2_human_influence",
  "cartnat_layer3_ecological_flow",
  "broad_leaf_type",
  "coniferous_leaf_type",
  "elevation",
  "forest_high_pressures_sum",
  "forest_low_pressures_sum",
  "forest_pressures_mode",
  "grassland",
  "road_fragmentation_average",
  "road_fragmentation_max",
  "mesh_fragmentation_average",
  "mesh_fragmentation_max",
  "small_woody_feature",
  "tree_cover_density",
  "permanent_water",
  "temporary_water",
  "any_water"
)

land.keep2 <- c(
  "id_carre_num",
  "broad_leaf_type",
  "coniferous_leaf_type",
  "elevation",
  "imperviousness",
  "grassland",
  "small_woody_feature",
  "tree_cover_density",
  "permanent_water",
  "temporary_water",
  "any_water"
)

land.stoc <-
  unique(stoc[, land.keep2])

# PCA on landscape features

pc.land.stoc <- dudi.pca(land.stoc[,-1],scannf = F, nf = 2)

# check eigenvalues 
screeplot(pc.land.stoc)

# test for number of axes
  #pca.land.testdim <- testdim(pc.land.stoc)
  #pca.land.testdim$nb
  #pca.land.testdim$nb.cor

s.corcircle(pc.land.stoc$co)

# PC1 : wetlands to drylands
# PC2 : urbanized / fragmented areas to forests

pc1.land.stoc <- cbind(stoc.id[, c("X", "Y")], pc.land.stoc$li[, 1])
colnames(pc1.land.stoc) = c("X", "Y", "PC1")
mapview(
  pc1.land.stoc,
  xcol = "X",
  ycol = "Y",
  zcol = "PC1",
  crs = 'epsg:27572'
)

pc2.land.stoc <- cbind(stoc.id[, c("X", "Y")], pc.land.stoc$li[, 2])
colnames(pc2.land.stoc) = c("X", "Y", "PC2")
mapview(
  pc2.land.stoc,
  xcol = "X",
  ycol = "Y",
  zcol = "PC2",
  crs = 'epsg:27572'
)

#---------------#
#### species ####
#---------------#

stoc.yr0 <- read.csv("data/stoc_yearly_abundances.csv")

# remove passer italiae (no trait data)
stoc.yr <- subset(stoc.yr0,birdlife_sci_name!="Passer italiae")

# add missing absences
rg.yr <- unique(stoc.yr$year)
carre.list <- unique(stoc.yr$id_carre)
species.list <- unique(stoc.yr$birdlife_sci_name)

stoc.index <- stoc.yr[,c("id_carre","year","birdlife_sci_name")]

stoc.index.all <- stoc.index %>% 
  count(id_carre, year,birdlife_sci_name) %>% 
  complete(id_carre, year,birdlife_sci_name, fill = list(n = 0))

stoc.all <- merge(stoc.yr,stoc.index.all,by = c("id_carre","year","birdlife_sci_name"),all = T)
stoc.all2 <- stoc.all[,c("id_carre","year","birdlife_sci_name","abundance")]
stoc.all2[is.na(stoc.all2$abundance),"abundance"] <- 0

# average abundance over years 
stoc.avg <-
  aggregate(
    stoc.all2$abundance,
    by = list(stoc.all2$id_carre, stoc.all2$birdlife_sci_name),
    FUN = "mean"
  )

colnames(stoc.avg) <- c("carre","species","avg.count")
stoc.avg$occurrence <- 0
stoc.avg[which(stoc.avg$avg.count > 0),"occurrence"] <- 1

# species ranking
occurrence.ranking <- aggregate(stoc.avg$occurrence,by = list(stoc.avg$species), FUN ="sum")
colnames(occurrence.ranking) <- c("species","total.count.stoc")

ggplot(occurrence.ranking) +
  geom_bar(aes(x = reorder(species, -total.count.stoc), y = total.count.stoc), stat = "identity")

occurrence.ranking <- occurrence.ranking[rev(order(occurrence.ranking[,2])),]
occurrence.ranking[1:50,]

