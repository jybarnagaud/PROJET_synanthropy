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
library(dplyr)
library(viridis)
library(PCAtest)
library(ape)
library(picante)
library(ggtree)
library(phylobase)
library(phylosignal)
library(mgcv)
library(MASS)
library(mgcViz)
library(ggeffects)
library(patchwork)
library(phylolm)

#------------#
#### data ####
#------------#

# acoustic traits

acoutr <-
  read.csv2("acoustic_traits_for_analysis.csv", row.names = 1)

# all other traits

stoc.traits <-
  read.csv2("outputs/full_trait_table.csv", row.names = 1)

# some grouping in categorical variables

stoc.traits0 <- stoc.traits
stoc.traits0[which(stoc.traits0$NomS == "Cinclus cinclus"), "Trophic.Niche"] <-
  "Invertivore"
stoc.traits0[which(stoc.traits0$NomS == "Alcedo atthis"), "Trophic.Niche"] <-
  "Omnivore"
stoc.traits0[which(stoc.traits0$NomS == "Lanius excubitor"), "Trophic.Niche"] <-
  "Invertivore"
stoc.traits0[which(stoc.traits0$X2019_uicn_redlist %in% c("NT", "VU")), "X2019_uicn_redlist"] <-
  "NT_VU"
stoc.traits0[which(stoc.traits0$X2019_uicn_redlist %in% c("Riverine", "Wetland")), "Habitat"] <-
  "River_Wetland"

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
    "syn1",
    "syn2",
    "syn3",
    "syn4"
  )]

#-------------------------------#
#### compute acoustic rarity ####
#-------------------------------#

# computing an euclidean distance matrix for acoustic traits

acou.dist <- as.matrix(dist(acoutr))

# false sites x species matrix for funrar

stsp.sim <- matrix(1, 3, nrow(acoutr))
colnames(stsp.sim) <- rownames(acoutr)

# functional rarity indices (non spatial)

acou.rar <- funrar(stsp.sim, acou.dist, rel_abund = F)

df.acou.rar <- data.frame(Ui = as.vector(acou.rar$Ui[, 2]),
                          Di = acou.rar$Di[1, ])

# all traits together

df.traits <-
  merge(
    eco.traits,
    df.acou.rar,
    by.x = "NomS",
    by.y = 0,
    all = T
  )

#----------------------------------#
#### multivariate "trait" space ####
#----------------------------------#

# gap filling with mean value

df.traits.nona <- df.traits %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))


# select traits to include

ts.data <-
  df.traits.nona[, !names(df.traits.nona) %in% c("NomS", "Affinity", "syn2", "syn3", "syn4")]
rownames(ts.data) <- df.traits.nona$NomS

# perform a PCA

pca.ts <- dudi.pca(ts.data, scannf = F, nf = 2)

# test for number of axes

pca.testdim <- testdim(pca.ts)
pca.testdim$nb
pca.testdim$nb.cor

# alternative test (for trial - similar result as testdim with Bonferroni procedure : 2 axes kept)
pcatest <- PCAtest(ts.data)

# check PCA

p1 <- fviz_pca_var(pca.ts)
p2 <- fviz_pca_ind(pca.ts)
p1 + p2

# where is affinity?

s.value(
  pca.ts$li,
  log(df.traits.nona$Affinity),
  method = "color",
  col = viridis_pal(alpha = 0.5)(7),
  ppoints = list(cex = 0.3)
)

# bivariate relationships between Affinity and PC axes

pc.af <-
  data.frame(
    pca.ts$li,
    Affinity = df.traits.nona$Affinity,
    LAffinity = log(df.traits.nona$Affinity)
  )

pc1.af <- ggplot(pc.af) +
  aes(x = Axis1, y = Affinity) +
  geom_point()

pc2.af <- ggplot(pc.af) +
  aes(x = Axis2, y = Affinity) +
  geom_point()

pc1.laf <- ggplot(pc.af) +
  aes(x = Axis1, y = LAffinity) +
  geom_point()

pc2.laf <- ggplot(pc.af) +
  aes(x = Axis2, y = LAffinity) +
  geom_point()

(pc1.af + pc2.af) / (pc1.laf + pc2.laf)

# bivariate relationships between Affinity and Ui, Di

laf.ui <- ggplot(df.traits) +
  aes(x = Ui, y = log(Affinity)) +
  geom_point()

laf.di <- ggplot(df.traits) +
  aes(x = Di, y = log(Affinity)) +
  geom_point()

laf.ui + laf.di

#----------------------------#
#### phylogenetic signals ####
#----------------------------#

# Thuiller's phylogenetic tree

wptree <- read.tree("data/phylogeny_birds_Thuiller2011.tre")
sp.in.phylo <- wptree$tip.label
sp.in.phylo <- sub("_", " ", sp.in.phylo)
wptree2 <- wptree
wptree2$tip.label <- sp.in.phylo

# implement taxonomic changes in the tree (see explore_trait_data.R for further details - code copy-pasted)

taxch <- read.csv2("data/taxonomic_update_thuiller_2011.csv")
sp.in.phylo2 <-
  data.frame(species.thuiller = sp.in.phylo, index = 1)
sp.in.phylo.match <-
  merge(sp.in.phylo2,
        taxch,
        by.x = "species.thuiller",
        by.y = "Thuiller_name",
        all = T)
sp.in.phylo.match[is.na(sp.in.phylo.match$STOC_name), "STOC_name"] <-
  sp.in.phylo.match[is.na(sp.in.phylo.match$STOC_name), "species.thuiller"]
rownames(sp.in.phylo.match) <- sp.in.phylo.match$species.thuiller
sp.in.phylo.match <- sp.in.phylo.match[wptree2$tip.label, ]

wptree3 <- wptree2
wptree3$tip.label <- sp.in.phylo.match$STOC_name

# prune tree

rownames(df.traits.nona) <- df.traits.nona$NomS
dif.tre <- setdiff(wptree3$tip.label, rownames(df.traits.nona))
afi.tree <- drop.tip(wptree3, dif.tre)

# check tree

ggtree(afi.tree, layout = 'circular') + geom_tiplab(size = 2, aes(angle =
                                                                    angle))

# phylogenetic signal

phy.data <- df.traits.nona
phy.data <- phy.data[afi.tree$tip.label, ]
phy.data <- as.data.frame(phy.data)

phylo4.data <- phylo4d(as(afi.tree, "phylo4"), phy.data)
barplot.phylo4d(phylo4.data, trait = c("Affinity", "syn1", "Ui", "Di"))
dotplot.phylo4d(phylo4.data, trait = c("Affinity", "syn1", "Ui", "Di"))

traits.physig <-
  phylo4.data[, c("Affinity",
                  "syn1",
                  "Ui",
                  "Di",
                  "Migration",
                  "RANGE_BIRDLIFE",
                  "SGIo")]
physig <- phyloSignal(traits.physig)
physig

#-------------------------------------------------------------#
#### model affinity against ecological and acoustic traits ####
#-------------------------------------------------------------#

# correlations among ecological traits

eco.tr <-
  df.traits.nona[, c("SGIo",
                     "RANGE_BIRDLIFE",
                     "Migration",
                     "LMass",
                     "Centroid.Latitude",
                     "syn1")]
chart.Correlation(eco.tr)
pc.eco <- dudi.pca(eco.tr, scannf = F, nf = 2)
p.eco <- fviz_pca_var(pc.eco)

df.mod <- df.traits.nona
df.mod$pc.eco1 <- pc.eco$li[, 1]
df.mod$pc.eco2 <- pc.eco$li[, 2]
df.mod$LUi <- log(df.mod$Ui)
df.mod$LDi <- log(df.mod$Di)

# correlation between Di and Ui and traits

chart.Correlation(df.mod[, c(
  "SGIo",
  "RANGE_BIRDLIFE",
  "Migration",
  "LMass",
  "Centroid.Latitude",
  "syn1",
  "Ui",
  "Di"
)])
chart.Correlation(df.mod[, c("pc.eco1", "pc.eco2", "Ui", "Di")])
chart.Correlation(df.mod[, c("pc.eco1", "pc.eco2", "LUi", "LDi")])

# log transform affinity for normalization

df.mod$LAffinity <- log(df.mod$Affinity)

# do the model

m.af <-
  gam(LAffinity ~ s(pc.eco1) + s(pc.eco2) + s(LUi) + s(LDi), data = df.mod)

# check concurvity

c.check1 <- concurvity(m.af, full = T)
c.check2 <- concurvity(m.af, full = F)

# check residuals

gam.check(m.af, old.style = T)

rsd <- residuals(m.af)

par(mfrow = c(2, 3))
qq.gam(m.af, rep = 100)
plot(fitted(m.af), rsd)
plot(df.mod$pc.eco1, rsd)
plot(df.mod$pc.eco2, rsd)
plot(df.mod$LUi, rsd)
plot(df.mod$LDi, rsd)

# check basis

par(mfrow = c(2, 2))
plot(m.af,
     residuals = TRUE,
     pch = 19,
     cex = .3)

# show model

summary(m.af)
par(mfrow = c(2, 2))
plot(m.af, residuals = T, cex = 5)

p.eco1.g <- ggpredict(m.af, "pc.eco1") %>%
  plot(residuals = T)

p.eco2.g <- ggpredict(m.af, "pc.eco2") %>%
  plot(residuals = T)

p.LUi.g <- ggpredict(m.af, "LUi") %>%
  plot(residuals = T)

p.LDi.g <- ggpredict(m.af, "LDi") %>%
  plot(residuals = T)

(p.eco1.g + p.eco2.g) / (p.LUi.g + p.LDi.g)

# refit as a linear model

m.af.l <-
  lm(LAffinity ~ pc.eco1 + pc.eco2 + LUi + LDi, data = df.mod)

summary(m.af.l)

p.eco1.l <- ggeffect(m.af.l, terms = "pc.eco1") %>%
  plot(residuals = T)

p.eco2.l <- ggeffect(m.af.l, terms = "pc.eco2")  %>%
  plot(residuals = T)

p.LUi.l <- ggeffect(m.af.l, terms = "LUi")  %>%
  plot(residuals = T)

p.LDi.l <- ggeffect(m.af.l, terms = "LDi")  %>%
  plot(residuals = T)

(p.eco1.l + p.eco2.l) / (p.LUi.l + p.LDi.l)

#--------------------------#
#### phylogenetic model ####
#--------------------------#

# note : just "to see" - no big reason for it because there is no signal in "affinity"

rs.l <- residuals(m.af.l)
rs.lp <- rs.l[(phylo4.data@data$NomS)]
phylo4.data@data$resid.lm <- rs.lp
barplot.phylo4d(phylo4.data, trait = "resid.lm")

test.sig.lm <- phyloSignal(phylo4.data[, c("resid.lm")])
test.sig.lm

# try a phylolm

df.mod.phy <- df.mod[afi.tree$tip.label, ]
m.af.p.BM <-
  phylolm(
    LAffinity ~ pc.eco1 + pc.eco2 + LUi + LDi,
    data = df.mod.phy,
    phy = afi.tree,
    model = "BM"
  )
m.af.p.lam <-
  phylolm(
    LAffinity ~ pc.eco1 + pc.eco2 + LUi + LDi,
    data = df.mod.phy,
    phy = afi.tree,
    model = "lambda"
  )
m.af.p.kap <-
  phylolm(
    LAffinity ~ pc.eco1 + pc.eco2 + LUi + LDi,
    data = df.mod.phy,
    phy = afi.tree,
    model = "kappa"
  )

# best fit model (not a good practice but no idea of any evolutionary model to choose)

AIC(m.af.p.BM)
AIC(m.af.p.lam)
AIC(m.af.p.kap)

# check residuals

par(mfrow = c(2, 2))
plot(m.af.p.lam)
plot(predict(m.af.p.lam), residuals(m.af.p.lam))
abline(h = 0)

# parameters (note : BM is provided just as a comparison because it is quite bad at residuals and AIC)

summary(m.af.l)
summary(m.af.p.lam)
summary(m.af.p.kap)
summary(m.af.p.BM)
