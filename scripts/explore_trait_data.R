library(ade4)
library(ggplot2)
library(funrar)
library(ape)
library(picante)
library(phytools)
library(ggtree)
library(viridis)
library(tidyverse)
library(phylosignal)
library(adephylo)
library(readxl)
library(phylobase)
library(openxlsx)
library(PCAtest)

#---------------------------------------#
#### data sur les traits acoustiques ####
#---------------------------------------#

acoutr2 <- read.csv2("acoustic_traits_for_analysis.csv",row.names = 1)

#-------------------------------------------------------#
#### Principal Component analysis on acoustic traits ####
#-------------------------------------------------------#

pca.acou <- dudi.pca(acoutr2,scannf=F,nf=5)

# check eigenvalues 
screeplot(pca.acou)

# test for number of axes
pca.acou.testdim <- testdim(pca.acou)
pca.acou.testdim$nb
pca.acou.testdim$nb.cor

# alternative test (for trial - similar result as testdim with Bonferroni procedure : 2 axes kept)
pcatest <- PCAtest(acoutr2)

# below we check four axes but only the two first will be kept for inference

# correlation circle 
s.corcircle(pca.acou$co,xax = 1, yax = 2)
s.corcircle(pca.acou$co,xax = 1, yax = 3)
s.corcircle(pca.acou$co,xax = 3, yax = 4)

s.label(pca.acou$li,xax = 1, yax = 2)
s.label(pca.acou$li,xax = 1, yax = 3)
s.label(pca.acou$li,xax = 3, yax = 4)

# interpretation : 
# 1 = dominant frequency
# 2 = spectral composition (~bandwidth)
# 3 = rythm
# 4 = song duration

#----------------------------#
#### phylogenetic signal #####
#----------------------------#

# Thuiller's phylogenetic tree 

wptree <- read.tree("data/phylogeny_birds_Thuiller2011.tre")
sp.in.phylo <- wptree$tip.label
sp.in.phylo <- sub("_", " ", sp.in.phylo)
wptree2 <- wptree
wptree2$tip.label <- sp.in.phylo

# implement taxonomic changes in the tree
# NOTE : Lanius meridionalis is absent from Thuiller's phylogeny. It has been assigned
# to the phylogenetic position of Lanius minor, which is absent from the STOC EPS.

taxch <- read.csv2("data/taxonomic_update_thuiller_2011.csv")
sp.in.phylo2 <- data.frame(species.thuiller = sp.in.phylo, index = 1)
sp.in.phylo.match <- merge(sp.in.phylo2,taxch,by.x = "species.thuiller", by.y = "Thuiller_name",all=T)
sp.in.phylo.match[is.na(sp.in.phylo.match$STOC_name),"STOC_name"] <- sp.in.phylo.match[is.na(sp.in.phylo.match$STOC_name),"species.thuiller"]  
rownames(sp.in.phylo.match) <- sp.in.phylo.match$species.thuiller 
sp.in.phylo.match <- sp.in.phylo.match[wptree2$tip.label,]

wptree3 <- wptree2
wptree3$tip.label <- sp.in.phylo.match$STOC_name

# prune tree

dif.tre <- setdiff(wptree3$tip.label,rownames(acoutr2))
acou.tree <- drop.tip(wptree3,dif.tre)

# check tree

ggtree(acou.tree, layout='circular') + geom_tiplab(size=2, aes(angle=angle))

#-----------------------------#
#### PCA axes on phylogeny ####
#-----------------------------#

pc.coord <- as_tibble(pca.acou$li)
pc.coord <-
  as_tibble(data.frame(label = rownames(pca.acou$li), pca.acou$li))

# plot trait space axes on phylogeny

source("functions/plotpcaxis.R")

plot.pcaxis(1)
ggsave("outputs/acoustic_pc1.png")

plot.pcaxis(2)
ggsave("outputs/acoustic_pc2.png")

plot.pcaxis(3)
ggsave("outputs/acoustic_pc3.png")

plot.pcaxis(4)
ggsave("outputs/acoustic_pc4.png")

# phylogenetic signal on pca axes

pc.coord.ord <- pca.acou$li
pc.coord.ord <- pc.coord.ord[acou.tree$tip.label,]
pc.coord.ord <- as.data.frame(pc.coord.ord)

phylo4.acou <- phylo4d(as(acou.tree,"phylo4"),pc.coord.ord)
barplot.phylo4d(phylo4.acou)
dotplot.phylo4d(phylo4.acou)
signal.acou <- phyloSignal(phylo4.acou)
signal.acou

# phylogenetic signal on raw traits

pc.coord.raw <- acoutr2
pc.coord.raw <- pc.coord.raw[acou.tree$tip.label,]
pc.coord.raw <- as.data.frame(pc.coord.raw)

phylo4.acou.raw <- phylo4d(as(acou.tree,"phylo4"),pc.coord.raw)
barplot.phylo4d(phylo4.acou.raw)
dotplot.phylo4d(phylo4.acou.raw)
signal.acou.raw <- phyloSignal(phylo4.acou.raw)
signal.acou.raw

# phylogenetically controlled PCA
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13515

phy.dist <- proxTips(acou.tree)
phy.pca <- ppca(phylo4.acou.raw, scannf = F, nfposi = 2, nfnega = 2)

scatter(phy.pca)
screeplot(phy.pca)
plot(phy.pca)
s.corcircle(phy.pca$c1)

# axis 1 = ~ spectral properties
# axis 2 = ~ frequency peaks

#---------------------------------------------------#
#### check the synanthropy index on test species ####
#---------------------------------------------------#

# trait data (from prepare_trait_data.R)

stoc.traits <- read.csv2("outputs/full_trait_table.csv",row.names = 1)

# distribution of synanthropy

chart.Correlation(stoc.traits[,c("syn1","syn2","syn3","syn4")])

sort.syn <- stoc.traits[,c("NomS","syn1","syn2","syn3","syn4")]

p.all.syn1 <- ggplot(sort.syn) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn1), stat = "identity")

most.syn1 <- sort.syn[rev(order(sort.syn$syn1)),]
most.syn2 <- sort.syn[rev(order(sort.syn$syn2)),]
most.syn3 <- sort.syn[rev(order(sort.syn$syn3)),]
most.syn4 <- sort.syn[rev(order(sort.syn$syn4)),]

# test species sets

common.species <- c(
  "Fringilla coelebs",
  "Passer domesticus",
  "Apus apus",
  "Sitta europaea",
  "Columba livia",
  "Anthus trivialis",
  "Oriolus oriolus",
  "Cuculus canorus",
  "Curruca communis",
  "Picus viridis"
)
stoc.traits.common <- subset(stoc.traits, NomS%in%common.species)

uncommon.species <- c(
  "Lullula arborea",
  "Dryobates minor",
  "Columba oenas",
  "Phylloscopus sibilatrix",
  "Passer montanus",
  "Corvus frugilegus",
  "Emberiza schoeniclus",
  "Pyrrhula pyrrhula",
  "Riparia riparia",
  "Alcedo atthis")
stoc.traits.uncommon <- subset(stoc.traits, NomS%in%uncommon.species)

rare.species <- c(
  "Fringilla montifringilla",
  "Apus pallidus",
  "Petronia petronia",
  "Anthus campestris",
  "Cinclus cinclus",
  "Coracias garrulus",
  "Emberiza hortulana",
  "Hippolais icterina",
  "Monticola solitarius",
  "Luscinia svecica"
)
stoc.traits.rare <- subset(stoc.traits, NomS%in%rare.species)

# look at the synanthropy index for common species

p.common.syn1 <- ggplot(stoc.traits.common) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn1), stat = "identity")
p.common.syn2 <- ggplot(stoc.traits.common) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn2), stat = "identity")
p.common.syn3 <- ggplot(stoc.traits.common) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn3), stat = "identity")
p.common.syn4 <- ggplot(stoc.traits.common) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn4), stat = "identity")

(p.common.syn1 + p.common.syn2) / (p.common.syn3 + p.common.syn4)
chart.Correlation(stoc.traits.common[,c("syn1","syn2","syn3","syn4")])

# "" for uncommon species

p.uncommon.syn1 <- ggplot(stoc.traits.uncommon) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn1), stat = "identity")
p.uncommon.syn2 <- ggplot(stoc.traits.uncommon) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn2), stat = "identity")
p.uncommon.syn3 <- ggplot(stoc.traits.uncommon) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn3), stat = "identity")
p.uncommon.syn4 <- ggplot(stoc.traits.uncommon) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn4), stat = "identity")

(p.uncommon.syn1 + p.uncommon.syn2) / (p.uncommon.syn3 + p.uncommon.syn4)
chart.Correlation(stoc.traits.uncommon[,c("syn1","syn2","syn3","syn4")])

# "" for rare species

p.rare.syn1 <- ggplot(stoc.traits.rare) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn1), stat = "identity")
p.rare.syn2 <- ggplot(stoc.traits.rare) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn2), stat = "identity")
p.rare.syn3 <- ggplot(stoc.traits.rare) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn3), stat = "identity")
p.rare.syn4 <- ggplot(stoc.traits.rare) +
  geom_bar(aes(x = reorder(NomS, -syn1), y = syn4), stat = "identity")

(p.rare.syn1 + p.rare.syn2) / (p.rare.syn3 + p.rare.syn4)
chart.Correlation(stoc.traits.rare[,c("syn1","syn2","syn3","syn4")])

