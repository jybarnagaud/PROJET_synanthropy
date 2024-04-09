#----------------------------------------------------------#
####       Exploration of Blackburn's affinity index    ####
# this is essentially a sub-analysis of explore_trait_data.R 
#----------------------------------------------------------#

library(funrar)
library(ade4)
library(factoextra)
library(PerformanceAnalytics)
library(adegraphics)
library(patchwork)

#------------#
#### data ####
#------------#

# acoustic traits

acoutr <- read.csv2("acoustic_traits_for_analysis.csv", row.names = 1)

# all other traits

stoc.traits <-
  read.csv2("outputs/full_trait_table.csv", row.names = 1)

# some grouping in categorical variables

stoc.traits0 <- stoc.traits
stoc.traits0[which(stoc.traits0$NomS == "Cinclus cinclus"),"Trophic.Niche"] <- "Invertivore"
stoc.traits0[which(stoc.traits0$NomS == "Alcedo atthis"),"Trophic.Niche"] <- "Omnivore"
stoc.traits0[which(stoc.traits0$NomS == "Lanius excubitor"),"Trophic.Niche"] <- "Invertivore"
stoc.traits0[which(stoc.traits0$X2019_uicn_redlist %in%c("NT","VU")),"X2019_uicn_redlist"] <- "NT_VU"
stoc.traits0[which(stoc.traits0$X2019_uicn_redlist %in%c("Riverine","Wetland")),"Habitat"] <- "River_Wetland"

# regularize body mass (otherwise weird results when doing an ordination)

stoc.traits0$LMass <- log(stoc.traits0$Mass)

# keep only some traits - note : Hand-Wing index removed because triggers weird patterns in ordination (Cisticola with large birds...)
eco.traits <-
  stoc.traits0[, c(
    "NomS",
    "SGIo",
    "RANGE_BIRDLIFE",
    "Affinity",
    "Migration",
    "LMass",
    "Centroid.Latitude",
    "syn1","syn2","syn3","syn4"
  )]

#-------------------------------#
#### compute acoustic rarity ####
#-------------------------------#

# computing an euclidean distance matrix for acoustic traits

acou.dist <- as.matrix(dist(acoutr))

# false sites x species matrix for funrar

stsp.sim <- matrix(1,3,nrow(acoutr))
colnames(stsp.sim) <- rownames(acoutr)

# functional rarity indices (non spatial)

acou.rar <- funrar(stsp.sim, acou.dist, rel_abund = F)

df.acou.rar <- data.frame(Ui = as.vector(acou.rar$Ui[,2]),
                          Di = acou.rar$Di[1,])

# all traits together

df.traits <- merge(eco.traits,df.acou.rar,by.x = "NomS",by.y = 0,all = T)

#----------------------------------#
#### multivariate "trait" space ####
#----------------------------------#

# gap filling with mean value

df.traits.nona <- df.traits %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) 


# select traits to include

ts.data <- df.traits.nona[,!names(df.traits.nona)%in%c("NomS","Affinity","syn2","syn3","syn4")]
rownames(ts.data) <- df.traits.nona$NomS

# perform a PCA

pca.ts <- dudi.pca(ts.data, scannf = F, nf = 2)

# test for number of axes

pca.testdim <- testdim(pca.ts)
pca.testdim$nb
pca.testdim$nb.cor

# alternative test (for trial - similar result as testdim with Bonferroni procedure : 2 axes kept)
pcatest <- PCAtest(ts.data.nona)

# check PCA

p1 <- fviz_pca_var(pca.ts)
p2 <- fviz_pca_ind(pca.ts)
p1 + p2

# where is affinity? 

s.value(pca.ts$li,log(df.traits.nona$Affinity),method = "color", col = viridis_pal(alpha = 0.5)(7),ppoints = list(cex = 0.3))

# bivariate relationships between Affinity and PC axes

pc.af <- data.frame(pca.ts$li,Affinity = df.traits.nona$Affinity,LAffinity = log(df.traits.nona$Affinity))

pc1.af <- ggplot(pc.af)+
  aes(x = Axis1, y= Affinity)+
  geom_point()

pc2.af <- ggplot(pc.af)+
  aes(x = Axis2, y= Affinity)+
  geom_point()

pc1.laf <- ggplot(pc.af)+
  aes(x = Axis1, y= LAffinity)+
  geom_point()

pc2.laf <- ggplot(pc.af)+
  aes(x = Axis2, y= LAffinity)+
  geom_point()

(pc1.af + pc2.af) / (pc1.laf + pc2.laf)

# bivariate relationships between Affinity and Ui, Di

laf.ui <- ggplot(df.traits)+
  aes(x = Ui, y = log(Affinity))+
  geom_point()

laf.di <- ggplot(df.traits)+
  aes(x = Di, y = log(Affinity))+
  geom_point()

laf.ui + laf.di
