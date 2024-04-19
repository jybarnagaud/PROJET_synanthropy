#---------------------------------------------------------------------------#
#### coinertia analysis between acoustic and non-acoustic traits #####
          # Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
                        # CESAB/ACOUCENE project
            # subproject on acoustic vs non acoustic traits
                              # april 2024
#---------------------------------------------------------------------------#

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

# functional - geographic - acoustic traits

stoc.traits <-
  as.data.frame(read_excel(
    "outputs/full_trait_table.xlsx",
    sheet = "all_traits",
    na = "NA"
  ))

# select traits 

traits0 <-
  stoc.traits[, !names(stoc.traits) %in% c(
    "nonphy.acou.PC1",
    "nonphy.acou.PC2",
    "phy.acou.PC1",
    "phy.acou.PC2"
  )]
rownames(traits0) <- traits0$NomS

# separate traits and covariates

traits <-
  traits0[, !names(traits0) %in% c("RANGE_BIRDLIFE", "Affinity", "X2019_uicn_redlist")]

# separate acoustic and functional traits

eco.traits <-
  traits[, c(
    "SGIo",
    "Hand-Wing.Index",
    "Beak.Depth",
    "Mass",
    "avg.clutch",
    "social",
    "longevity",
    "Migration",
    "Habitat",
    "Trophic.Niche",
    "Primary.Lifestyle",
    "Centroid.Latitude",
    "Range.Size",
    "syn1",
    "syn2",
    "syn3",
    "syn4"
  )]

acou.traits <-
  traits[, c("Ht","Hf","nfp","nfp.cqv","pf",
             "pf.cqv","actfract","sp.bdw","ugof.iso","npvi.ioi","npvi.ioi.cqv",
             "n.syll","dur","syll.dur.med","ioi.dur.med","syll.dur.med.cqv")]

covariates <-
  traits0[, c("RANGE_BIRDLIFE", "Affinity", "X2019_uicn_redlist")]

# separate functional and geographic traits

fun.traits <-
  traits[, c(
    "Hand-Wing.Index",
    "Beak.Depth",
    "Mass",
    "avg.clutch",
    "social",
    "longevity",
    "Migration",
    "Habitat",
    "Trophic.Niche",
    "Primary.Lifestyle"
  )]

geo.traits <-
  traits[, c("SGIo", "Centroid.Latitude", "Range.Size", "syn1", "syn2","syn3","syn4")]

# "conservation" traits

cons.traits <-
  traits0[, names(traits0) %in% c("RANGE_BIRDLIFE", "Affinity", "X2019_uicn_redlist")]

cons.traits$X2019_uicn_redlist <-
  factor(cons.traits$X2019_uicn_redlist)

# gap filling (mean) for missing values

geo.traits.filled <- geo.traits %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE))))

fun.traits.filled <- fun.traits %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE))))

cons.traits.filled <- cons.traits %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE))))

## correlation matrices --------------------------------------------------------

chart.Correlation(fun.traits[,c(
  "Hand-Wing.Index",
  "Beak.Depth",
  "Mass",
  "Migration",
  "avg.clutch",
  "social",
  "longevity")])
  
chart.Correlation(geo.traits)

## matrices for ordinations ----------------------------------------------------

# keep only syn4 due to high correlation among all "syn"

geo.traits.sub <- geo.traits.filled[,!names(geo.traits.filled)%in%c("syn1","syn2","syn3")]

# match filled trait matrices

trait.filled <- merge(fun.traits.filled,geo.traits.sub,by=0)

## single - table ordination on functional traits ------------------------------

trait.ord <- trait.filled[,-1]
rownames(trait.ord) <- trait.filled$Row.names

trait.ord$Habitat <- factor(trait.ord$Habitat)
trait.ord$Trophic.Niche <- factor(trait.ord$Trophic.Niche)
trait.ord$Primary.Lifestyle <- factor(trait.ord$Primary.Lifestyle)

hs.trait <- dudi.hillsmith(trait.ord,scannf = F, nf = 2)

# ordination plots

screeplot(hs.trait)
100* hs.trait$eig/sum(hs.trait$eig)
100* cumsum(hs.trait$eig)/sum(hs.trait$eig)

# no test of the number of axes available for hillsmith. done on a PCA

trait.quanti <- trait.ord

trait.quanti$Habitat <-
  as.numeric(factor(
    trait.quanti$Habitat,
    labels = c(
      "Rock",
      "Human Modified",
      "Grassland",
      "Shrubland",
      "Wetland",
      "Riverine",
      "Woodland",
      "Forest"
    )
  ))

trait.quanti$Trophic.Niche <-
  as.numeric(factor(
    trait.quanti$Trophic.Niche,
    levels =
      c(
        "Vertivore",
        "Aquatic predator",
        "Invertivore",
        "Omnivore",
        "Granivore"
      )
  ))


trait.quanti$Primary.Lifestyle <-
  as.numeric(factor(
    trait.quanti$Primary.Lifestyle,
    levels = c("Terrestrial", "Insessorial", "Generalist", "Aerial")
  ))

pc.trait.quanti <- dudi.pca(trait.quanti,scannf = F, nf = 2 )

# Dray's test for number of axes

pc.trait.quanti.test <- testdim(pc.trait.quanti)
pc.trait.quanti.test$nb
pc.trait.quanti.test$nb.cor

# look at the Hill & Smith analysis with two axes

s.label(
  hs.trait$li,
  plabels = list(
    cex = 0.6,
    col = "gray60",
    boxes = list(draw = F)
  ),
  ppoints = list(cex = 0),
  pgrid = list(draw = F)
)

s.arrow(
  hs.trait$co,
  add = T,
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  )

)

# zoom into variables

p1 <- s.arrow(
  hs.trait$co,
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  ),
  pgrid = list(draw = F)
)

p2 <- s.label(
  hs.trait$li,
  plabels = list(
    cex = 0.5,
    boxes = list(draw = F)
  ),
  ppoints = list(cex = 0),
  pgrid = list(draw = F)
)

# png("outputs/hillsmith-functional-traits.png", width = 3000, height = 1500, res = 300)

cbindADEg(p1, p2, plot = TRUE)

# dev.off()

## single - table ordination on acoustic traits --------------------------------

pca.acou <- dudi.pca(acou.traits,scannf=F,nf=4)

# ordination plots

screeplot(pca.acou)

pca.acou$eig/sum(pca.acou$eig)
cumsum(pca.acou$eig)/sum(pca.acou$eig)

# test for number of axes

pca.acou.testdim <- testdim(pca.acou)
pca.acou.testdim$nb
pca.acou.testdim$nb.cor

# pcatest <- PCAtest(acou.traits)


# ordination plots

p3.1 <- s.corcircle(
  pca.acou$co,
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  ),
  pgrid = list(draw = F)
)

p3.2 <- s.corcircle(
  pca.acou$co,
  xax = 2,
  yax = 3,
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  ),
  pgrid = list(draw = F)
)

p3.3 <- s.corcircle(
  pca.acou$co,
  xax = 3,
  yax = 4,
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  ),
  pgrid = list(draw = F)
)


p4.1 <- s.label(
  pca.acou$li,
  plabels = list(
    cex = 0.5,
    alpha = 0.5,
    col = "gray30",
    boxes = list(draw = F)
  ),
  ppoints = list(cex = 0),
  pgrid = list(draw = F)
)

p4.2 <- s.label(
  pca.acou$li,
  xax = 2, yax = 3,
  plabels = list(
    cex = 0.5,
    alpha = 0.5,
    col = "gray30",
    boxes = list(draw = F)
  ),
  ppoints = list(cex = 0),
  pgrid = list(draw = F)
)

p4.3 <- s.label(
  pca.acou$li,
  xax = 3, yax = 4,
  plabels = list(
    cex = 0.5,
    alpha = 0.5,
    col = "gray30",
    boxes = list(draw = F)
  ),
  ppoints = list(cex = 0),
  pgrid = list(draw = F)
)

#png("outputs/pca-acoustic-traits-variables.png", width = 3000, height = 1500, res = 300)

cbindADEg(p3.1,p3.2,p3.3, plot = TRUE)

#dev.off()

#png("outputs/pca-acoustic-traits-species.png", width = 3000, height = 1500, res = 300)

cbindADEg(p4.1,p4.2,p4.3, plot = TRUE)

#dev.off()

## coinertia analysis ----------------------------------------------------------

coi.acou.fun <- coinertia(hs.trait, pca.acou, scannf = F, nf = 2)

# permutation test (with the quantitative matrix - not available for Hill Smith)

coi.acou.fun.quant <- coinertia(pc.trait.quanti, pca.acou, scannf = F, nf = 2)
coi.test <- randtest(coi.acou.fun.quant)
plot(coi.test)

# ordination plots

screeplot(coi.acou.fun)
coi.acou.fun$eig/sum(coi.acou.fun$eig)
cumsum(coi.acou.fun$eig)/sum(coi.acou.fun$eig)

p5 <- s.label(
  coi.acou.fun$co ,
  ppoints = list(cex= 0),
  plabels = list(
    cex = 0.7,
    col = "darkred",
    boxes = list(draw = F)
  ),
  pgrid = list(draw = F)
)

p6 <- s.label(
  coi.acou.fun$li,
  ppoints = list(cex= 0),
  plabels = list(
    cex = 0.7,
    col = "darkblue",
    boxes = list(draw = F)
  )
)

p56=superpose(p5,p6)


p7 <- plot(
  coi.acou.fun,
  plabels = list(
    cex = 0.5,
    alpha = 0.5,
    col = "gray30",
    boxes = list(draw = F)),
  pgrid = list(draw = F))[[4]]

png("outputs/coinertia-functional-acoustic-traits.png", width = 3000, height = 1500, res = 300)

cbindADEg(p56, p7, plot = TRUE)

dev.off()


##  acoustic - conservation variables ------------------------------------------

df.acou <- data.frame(pca.acou$li,
                      cons.traits.filled)
colnames(df.acou) <- c("acou.pc1","acou.pc2","acou.pc3","acou.pc4","rangesize","attract","redlist")

df.acou$lattract <- log(df.acou$attract)
df.acou$lrangesize <- log(df.acou$rangesize)

p11 <- ggplot(df.acou)+
  aes(x = acou.pc1,y=lattract)+
  geom_point()

p12 <- ggplot(df.acou)+
  aes(x = acou.pc2,y=lattract)+
  geom_point()

p13 <- ggplot(df.acou)+
  aes(x = acou.pc3,y=lattract)+
  geom_point()

p14 <- ggplot(df.acou)+
  aes(x = acou.pc4,y=lattract)+
  geom_point()

(p11 + p12) / (p13 + p14) 

p15 <- ggplot(df.acou)+
  aes(x = acou.pc1,y=lrangesize)+
  geom_point()

p16 <- ggplot(df.acou)+
  aes(x = acou.pc2,y=lrangesize)+
  geom_point()

p17 <- ggplot(df.acou)+
  aes(x = acou.pc3,y=lrangesize)+
  geom_point()

p18 <- ggplot(df.acou)+
  aes(x = acou.pc4,y=lrangesize)+
  geom_point()

(p15 + p16) / (p17 + p18) 

