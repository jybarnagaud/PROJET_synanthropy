#-----------------------------------------------------#
#### analysis of acoustic uniqueness  #####
# Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
              # CESAB/ACOUCENE project
      # subproject on acoustic vs non acoustic traits
                  # april 2024
#-----------------------------------------------------#

library(readxl)
library(dplyr)
library(tidyr)
library(viridis)
library(ade4)
library(adegraphics)
library(factoextra)
library(PCAtest)
library(PerformanceAnalytics)
library(patchwork)

## data ------------------------------------------------------------------------
    
stoc.traits <-
  as.data.frame(read_excel(
    "outputs/full_trait_table.xlsx",
    sheet = "all_traits",
    na = "NA"
  ))

acou.traits <-
  traits[, c("Ht","Hf","nfp","nfp.cqv","pf",
             "pf.cqv","actfract","sp.bdw","ugof.iso","npvi.ioi","npvi.ioi.cqv",
             "n.syll","dur","syll.dur.med","ioi.dur.med","syll.dur.med.cqv")]

## functional rarity -----------------------------------------------------------
    
# computing an euclidean distance matrix for acoustic traits
    
acou.dist <- as.matrix(dist(acou.traits))
    
# false sites x species matrix for funrar
    
stsp.sim <- matrix(1,3,nrow(acou.traits))
colnames(stsp.sim) <- rownames(acou.traits)
    
# functional rarity indices (non spatial)
    
acou.rar <- funrar(stsp.sim, acou.dist, rel_abund = F)
    
df.acou.rar <- data.frame(Ui = as.vector(acou.rar$Ui[,2]),
                              Di = acou.rar$Di[1,])


# acoustic PCA with rarity

pca.acou <- dudi.pca(acou.traits,scannf=F,nf=4)

s.value(pca.acou$li,z = df.acou.rar$Ui,method = "color")

df.acou.rar$acou.pc1 <- pca.acou$li[,1]
df.acou.rar$acou.pc2 <- pca.acou$li[,2]
df.acou.rar$acou.pc3 <- pca.acou$li[,3]
df.acou.rar$acou.pc4 <- pca.acou$li[,4]

p1 <- ggplot(df.acou.rar)+
  aes(x = acou.pc1,y=Ui)+
  geom_point()

p2 <- ggplot(df.acou.rar)+
  aes(x = acou.pc2,y=Ui)+
  geom_point()

p3 <- ggplot(df.acou.rar)+
  aes(x = acou.pc3,y=Ui)+
  geom_point()

p4 <- ggplot(df.acou.rar)+
  aes(x = acou.pc4,y=Ui)+
  geom_point()

(p1 + p2) / (p3 + p4)

ggsave("outputs/acoustic_uniqueness_vs_pca.png")


p5 <- ggplot(df.acou.rar)+
  aes(x = acou.pc1,y=Di)+
  geom_point()

p6 <- ggplot(df.acou.rar)+
  aes(x = acou.pc2,y=Di)+
  geom_point()

p7 <- ggplot(df.acou.rar)+
  aes(x = acou.pc3,y=Di)+
  geom_point()

p8 <- ggplot(df.acou.rar)+
  aes(x = acou.pc4,y=Di)+
  geom_point()

(p5 + p6) / (p7 + p8)

ggsave("outputs/acoustic_distinctiveness_vs_pca.png")


