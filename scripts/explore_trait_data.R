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

setwd("C:/Users/jeany/OneDrive/Documents/PROJET_synanthropy/data")

#----------------------------------------------------------------#
#### acoustic traits retrieved by M. Busana from Jean Roché   ####
#----------------------------------------------------------------#

acoutr <- read.csv("cleaned_maad_4_histo.csv")
rownames(acoutr) <- acoutr$birdlife_sci_name

# gap filling : centroïd values for two NA on "duration songs" (species for which this trait is irrelevant / not computable)
acoutr[is.na(acoutr$duration_song), "duration_song"] <-
  mean(acoutr$duration_song, na.rm = T)

# keep
acoutr2 <-
  acoutr[, c(
    "duration_song",
    "centroid_f",
    "Ht_temporal_entropy",
    "LFC_spectral_cover",
    "MFC_spectral_cover",
    "HFC_spectral_cover",
    "Hf_spectral_entropy",
    "number_of_freq_peaks",
    "Hf_Havrda_spectral_entropy",
    "Hf_Renyi_spectral_entropy",
    "Hf_pairedShannon_spectral_entropy",
    "Hf_gamma_spectral_entropy",
    "Hf_GiniSimpson_spectral_entropy",
    "peak_freq",#
    "peak_freq_amp_by_mean_amp",
    "peak_f_roi",
    "centroid_f_roi",
    "num_syllables_per_unit_time",
    "syllable_duration_mean",
    "syllable_duration_std",
    "mean_gap_duration",
    "gap_duration_std",
    "sound_per_silence_ratio",
    "sound_per_vocalization_ratio"
  )]


acoutr2 <-
  acoutr[, c(
    "duration_song",
    "centroid_f",
    "LFC_spectral_cover",
    "MFC_spectral_cover",
    "HFC_spectral_cover",
    "Hf_spectral_entropy",
    "number_of_freq_peaks",
    "peak_freq",
    "syllable_duration_mean",
    "sound_per_vocalization_ratio"
  )]

#-------------------------------------------------------#
#### Principal Component analysis on acoustic traits ####
#-------------------------------------------------------#

pca.acou <- dudi.pca(acoutr2,scannf=F,nf=5)

# correlation circle 
s.corcircle(pca.acou$co,xax = 1, yax = 2)
s.corcircle(pca.acou$co,xax = 1, yax = 3)
s.corcircle(pca.acou$co,xax = 3, yax = 4)

s.label(pca.acou$li,xax = 1, yax = 2)
s.label(pca.acou$li,xax = 1, yax = 3)
s.label(pca.acou$li,xax = 3, yax = 4)


# 4 = durée totale du chant
# 3 = rythme
# 2 = composition spectrale (~bandwidth)
# 1 = fréquence dominante

#----------------------------#
#### phylogenetic signal #####
#----------------------------#

# Thuiller's phylogenetic tree 

wptree <- read.tree("phylogeny_birds_Thuiller2011.tre")
sp.in.phylo <- wptree$tip.label
sp.in.phylo <- sub("_", " ", sp.in.phylo)
wptree2 <- wptree
wptree2$tip.label <- sp.in.phylo

# implement taxonomic changes in the tree
# NOTE : Lanius meridionalis is absent from Thuiller's phylogeny. It has been assigned
# to the phylogenetic position of Lanius minor, which is absent from the STOC EPS.

taxch <- read.csv2("taxonomic_update_thuiller_2011.csv")
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
ggtree(acou.tree, layout='circular') + geom_tiplab(size=2, aes(angle=angle))

#-----------------------------#
#### PCA axes on phylogeny ####
#-----------------------------#

pc.coord <- as_tibble(pca.acou$li)
pc.coord <-
  as_tibble(data.frame(label = rownames(pca.acou$li), pca.acou$li))

# plot trait space axes on phylogeny

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

plot.pcaxis(1)
ggsave("acoustic_pc1.png")

plot.pcaxis(2)
ggsave("acoustic_pc2.png")

plot.pcaxis(3)
ggsave("acoustic_pc3.png")

plot.pcaxis(4)
ggsave("acoustic_pc4.png")

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
phy.pca <- ppca(phylo4.acou.raw)

scatter(phy.pca)
screeplot(phy.pca)
plot(phy.pca)

#-----------------------------------------------#
##### non acoustic traits : hand-wing index #####
#-----------------------------------------------#

avonet <-
  as.data.frame(read_excel("AVONET Supplementary dataset 1.xlsx", sheet = "AVONET1_BirdLife"))

hwi <- avonet[, c("Species1", "Hand-Wing.Index")]
hwi.sub <- subset(hwi, Species1 %in% acou.tree$tip.label)

missp.hwi <- setdiff(acou.tree$tip.label, hwi$Species1)

missp.hwi.thuiller <-
  c(
    "Curruca conspicillata",
    "Curruca undata",
    "Curruca cantillans",
    "Curruca melanocephala",
    "Curruca communis",
    "Curruca curruca",
    "Curruca hortensis",
    "Luscinia svecica"
  )

missp.hwi.avonet <-
  c(
    "Sylvia conspicillata",
    "Sylvia undata",
    "Sylvia cantillans",
    "Sylvia melanocephala",
    "Sylvia communis",
    "Sylvia curruca",
    "Sylvia hortensis",
    "Cyanecula svecica"
  )

tax.hwi <- data.frame(missp.hwi.thuiller, missp.hwi.avonet)


hwi2 <-
  merge(hwi,
        tax.hwi,
        by.x = "Species1",
        by.y = "missp.hwi.avonet",
        all = T)
hwi2[!is.na(hwi2$missp.hwi.thuiller), "Species1"] <-
  hwi2[!is.na(hwi2$missp.hwi.thuiller), "missp.hwi.thuiller"]

hwi.sub <- subset(hwi2, Species1 %in% acou.tree$tip.label)[,-3]

#----------------------------------------------#
##### non acoustic traits : affinity index #####
#----------------------------------------------#

afi <- as.data.frame(read_excel(
  "xenocanto_acoucene.xlsx",
  sheet = "xeno canto europe au 25fev2020"
))

afi.sub <- subset(afi, Nom_scientifique %in% acou.tree$tip.label)

missp.afi <- setdiff(acou.tree$tip.label, afi$Nom_scientifique)

missp.afi.thuiller <-
  c(
    "Curruca conspicillata",
    "Curruca undata",
    "Curruca cantillans",
    "Curruca melanocephala",
    "Curruca communis",
    "Curruca curruca",
    "Curruca hortensis",
    "Leiopicus medius",
    "Corvus monedula",
    "Saxicola torquatus" 
  )

missp.afi.source <-
  c(
    "Sylvia conspicillata",
    "Sylvia undata",
    "Sylvia cantillans",
    "Sylvia melanocephala",
    "Sylvia communis",
    "Sylvia curruca",
    "Sylvia hortensis",
    "Dendrocopos medius",
    "Coloeus monedula",
    "Saxicola torquata"
  )

tax.afi <- data.frame(missp.afi.thuiller, missp.afi.source)

afi2 <-
  merge(afi,
        tax.afi,
        by.x = "Nom_scientifique",
        by.y = "missp.afi.source",
        all = T)
afi2[!is.na(afi2$missp.afi.thuiller), "Nom_scientifique"] <-
  afi2[!is.na(afi2$missp.afi.thuiller), "missp.afi.thuiller"]

afi.sub <- subset(afi2, Nom_scientifique %in% acou.tree$tip.label)[,c("Nom_scientifique","No")]
colnames(afi.sub) <- c("Species1","Affinity")

#---------------------------------------------#
##### non acoustic traits : habitat index #####
#---------------------------------------------#

hab <- as.data.frame(read_excel(
  "SGI Godet et al.xlsx",
  sheet = "SXI-FINALv4"
))

taxstoc <- as.data.frame(read_excel(
  "referentiel stoceps.xlsx",
  sheet = "TEEspèces"
))
  
taxupdate <- read.csv2(
  "reftax_STOCEPS_thuiller.csv"
)

taxstoc.update <- merge(taxstoc,taxupdate,by.x = "Espèce",by.y = "code_STOCEPS",all=T)
missp.habtax <- setdiff(acou.tree$tip.label, hab2$NomS)

habtax <-
  merge(
    hab,
    taxstoc.update,
    by.y = "Espèce",
    by.x = "ESPECE",
    all.y = F,
    all.x = T
  )
habtax[!is.na(habtax$species_Thuiller), "NomS"] <-
  habtax[!is.na(habtax$species_Thuiller), "species_Thuiller"]

hab2 <- subset(habtax,NomS %in% acou.tree$tip.label)[,c("NomS","SGIo","RANGE_BIRDLIFE")]

#------------------------------------------------#
#### non acoustic traits : UICN threat status ####
#------------------------------------------------#

iucn <- read.csv2("uicn_threat_status_2019.csv")
iucn.hab <- 
missp.iucn <- setdiff(acou.tree$tip.label, iucn$species) # same missing species as for hwi

iucn2 <-
  merge(iucn,
        tax.hwi,
        by.x = "species",
        by.y = "missp.hwi.avonet",
        all = T)
iucn2[!is.na(iucn2$missp.hwi.thuiller), "species"] <-
  iucn2[!is.na(iucn2$missp.hwi.thuiller), "missp.hwi.thuiller"]

iucn.sub <- subset(iucn2, species %in% acou.tree$tip.label)[,-3]

#-------------------------------------------------------#
#### match all non acoustic traits in a single table ####
#-------------------------------------------------------#

habi.afi <-
  merge(hab2,
        afi.sub,
        by.x = "NomS",
        by.y = "Species1",
        all = T)

habi.afi.hwi <-
  merge(habi.afi,
        hwi.sub,
        by.x = "NomS",
        by.y  = "Species1",
        all = T)

habi.afi.hwi.iucn <-
  merge(habi.afi.hwi,
        iucn.sub,
        by.x = "NomS",
        by.y  = "species",
        all = T)

non.acou.traits 