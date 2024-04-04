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
