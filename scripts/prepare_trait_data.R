#-------------------------------------------------------------#
#### Compile and prepare acoustic and non-acoustic traits #####
    # Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
                  # CESAB/ACOUCENE project
        # subproject on acoustic vs non acoustic traits
                        # april 2024
#-------------------------------------------------------------#

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

#----------------------------------------------------------------#
#### acoustic traits retrieved by M. Busana from Jean Roché   ####
#----------------------------------------------------------------#

acoutr0 <-read.csv(
    "data/maad_features_single_row_by_species_and_no_drum.csv",
  )

traits.keep <- c(
  "Ht_temporal_entropy",
  "Hf_spectral_entropy",
  "number_of_freq_peaks",
  "number_of_freq_peaks_cqv",
  "peak_freq",
  "peak_freq_cqv",
  "ACTfract_mean_proportion_of_points_above_6dB_threshold",
  "spectral_bandwidth_90",
  "ugof_isochronous",
  "npvi_ioi",
  "npvi_ioi_cqv",
  "num_syllables_per_unit_time",
  "duration_song",
  "syllable_duration_median",
  "ioi_duration_median",
  "syllable_duration_median_cqv"
)

# duplicated species : average over traits of interest
acoutr1 <-
  aggregate(acoutr0[, traits.keep],
            by = list(acoutr0$birdlife_sci_name),
            FUN = "mean")

# simpler acronyms

colnames(acoutr1)[-1] <- c("Ht","Hf","nfp","nfp.cqv","pf",
                           "pf.cqv","actfract","sp.bdw","ugof.iso","npvi.ioi","npvi.ioi.cqv",
                           "n.syll","dur","syll.dur.med","ioi.dur.med","syll.dur.med.cqv")

# species as rownames (for ade4)
rownames(acoutr1) <- acoutr1$Group.1
acoutr <- acoutr1[,-1]

# gap filling : centroïd values for two NA on "duration songs" (species for which this trait is irrelevant / not computable)

acoutr <- acoutr %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE))))

# save trait data
write.csv2(acoutr,"acoustic_traits_for_analysis.csv")

#-------------------------------------------------------#
#### Principal Component analysis on acoustic traits ####
#-------------------------------------------------------#

pca.acou <- dudi.pca(acoutr,scannf=F,nf=5)

screeplot(pca.acou)

# we keep the two first PC axes - see explore_trait_data.R for choice material

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

dif.tre <- setdiff(wptree3$tip.label,rownames(acoutr))
acou.tree <- drop.tip(wptree3,dif.tre)

#-----------------------------#
#### PCA axes on phylogeny ####
#-----------------------------#

pc.coord <- as_tibble(pca.acou$li)
pc.coord <-
  as_tibble(data.frame(label = rownames(pca.acou$li), pca.acou$li))

# phylogenetic signal on pca axes

pc.coord.ord <- pca.acou$li
pc.coord.ord <- pc.coord.ord[acou.tree$tip.label,]
pc.coord.ord <- as.data.frame(pc.coord.ord)

phylo4.acou <- phylo4d(as(acou.tree,"phylo4"),pc.coord.ord)

# phylogenetic signal on raw traits

pc.coord.raw <- acoutr
pc.coord.raw <- pc.coord.raw[acou.tree$tip.label,]
pc.coord.raw <- as.data.frame(pc.coord.raw)

phylo4.acou.raw <- phylo4d(as(acou.tree,"phylo4"),pc.coord.raw)

# phylogenetically controlled PCA
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13515

phy.dist <- proxTips(acou.tree)
phy.pca <- ppca(phylo4.acou.raw, scannf = F, nfposi = 2, nfnega = 2)

# axis 1 = ~ spectral properties
# axis 2 = ~ frequency peaks

#---------------------------------------------------------------#
##### non acoustic traits : "functional traits" from AVONET #####
#---------------------------------------------------------------#

# note : also includes range size and latitudinal centroid (geographic traits)
avonet <-
  as.data.frame(read_excel("data/AVONET Supplementary dataset 1.xlsx", sheet = "AVONET1_BirdLife",na="NA"))

funtr <- avonet[, c("Species1", "Hand-Wing.Index","Beak.Depth", "Mass","Migration","Habitat","Trophic.Niche","Primary.Lifestyle","Centroid.Latitude", "Range.Size")]
funtr.sub <- subset(funtr, Species1 %in% acou.tree$tip.label)

# taxonomic homogenization

missp.avonet <- setdiff(acou.tree$tip.label, funtr$Species1)

missp.avonet.thuiller <-
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

missp.avonet.avonet <-
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

tax.funtr <- data.frame(missp.avonet.thuiller, missp.avonet.avonet)

# functional traits

funtr2 <-
  merge(funtr,
        tax.funtr,
        by.x = "Species1",
        by.y = "missp.avonet.avonet",
        all = T)
funtr2[!is.na(funtr2$missp.avonet.thuiller), "Species1"] <-
  funtr2[!is.na(funtr2$missp.avonet.thuiller), "missp.avonet.thuiller"]

funtr.sub0 <- subset(funtr2, Species1 %in% acou.tree$tip.label)
funtr.sub1 <- funtr.sub0[,!names(funtr.sub0) == "missp.avonet.thuiller"]

#------------------------------------------------------------------------------#
##### non acoustic traits : repro and social from C. Sekercioglu's database ####
#------------------------------------------------------------------------------#

# !! use and distribution of these data are subject to Cagan Sekercioglu's agreement (cagan1@gmail.com)
# !! non homogeneous data and/or decimal separation in this file

breed0 <- read.csv2("data/2021 Barnagaud Birdbase TRAVAIL 2022.csv",sep=";",dec=".")
breed <- breed0[,c("Gill.et.al.2020..IOC.World.Bird.List..v10.2.","Clutch.Ave","Social","Known.Longevity")]
colnames(breed) <- c("Species1","avg.clutch","social","longevity")

# match taxonomy

breed["Dendrocoptes medius","Species1"] <- "Leiopicus medius"
breed["Coloeus monedula","Species1"] <- "Corvus monedula"

# subset species

breed1 <- subset(breed,Species1%in%funtr.sub1$Species1)

# rework trait "social" : average of all possible classes

soc.breed <- breed1$social
soc.breed.split <- strsplit(soc.breed,",")

soc.breed.av <- unlist(lapply(soc.breed.split,FUN=function(x){mean(as.numeric(x))}))

breed2 <- breed1
breed2$social <- soc.breed.av

# change longevity to numeric

breed2$longevity <- as.numeric(breed2$longevity)

# all functional traits

funtr.sub <- merge(funtr.sub1,breed2,by="Species1")

#----------------------------------------------#
##### non acoustic traits : affinity index #####
#----------------------------------------------#

afi <- as.data.frame(read_excel(
  "data/xenocanto_acoucene.xlsx",
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
  "data/SGI Godet et al.xlsx",
  sheet = "SXI-FINALv4"
))

taxstoc <- as.data.frame(read_excel(
  "data/referentiel stoceps.xlsx",
  sheet = "TEEspèces"
))
  
taxupdate <- read.csv2(
  "data/reftax_STOCEPS_thuiller.csv"
)

taxstoc.update <- merge(taxstoc,taxupdate,by.x = "Espèce",by.y = "code_STOCEPS",all=T)


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

missp.habtax <- setdiff(acou.tree$tip.label, hab2$NomS)

#------------------------------------------------#
#### non acoustic traits : UICN threat status ####
#------------------------------------------------#

iucn <- read.csv2("data/uicn_threat_status_2019.csv")
iucn.hab <- 
missp.iucn <- setdiff(acou.tree$tip.label, iucn$species) # same missing species as for functional traits

iucn2 <-
  merge(iucn,
        tax.funtr,
        by.x = "species",
        by.y = "missp.avonet.avonet",
        all = T)
iucn2[!is.na(iucn2$missp.avonet.thuiller), "species"] <-
  iucn2[!is.na(iucn2$missp.avonet.thuiller), "missp.avonet.thuiller"]

iucn.sub <- subset(iucn2, species %in% acou.tree$tip.label)[,-3]

#-----------------------------------------------#
#### non acoustic traits : synanthropy index ####
#-----------------------------------------------#

# synanthropy index computed by Nicolas Casajus (march 2024) from CARTNAT

syn <- read.csv("data/acoucene_ssi_results_121species.csv")

# separate synanthropy indices from different CARTNAT layers

syn.cartnat1 <- subset(syn,cartnat_layer == 1)
syn.cartnat2 <- subset(syn,cartnat_layer == 2)
syn.cartnat3 <- subset(syn,cartnat_layer == 3)
syn.cartnat4 <- subset(syn,cartnat_layer == 4)

syn.wide <- as.data.frame(pivot_wider(syn,names_from = cartnat_layer, values_from = ssi_mean))
colnames(syn.wide) <- c("species", "syn1","syn2","syn3","syn4")

#-------------------------------------------------------#
#### match all non acoustic traits in a single table ####
#-------------------------------------------------------#

habi.afi <-
  merge(hab2,
        afi.sub,
        by.x = "NomS",
        by.y = "Species1",
        all = T)

habi.afi.funtr <-
  merge(habi.afi,
        funtr.sub,
        by.x = "NomS",
        by.y  = "Species1",
        all = T)

habi.afi.funtr.iucn <-
  merge(habi.afi.funtr,
        iucn.sub,
        by.x = "NomS",
        by.y  = "species",
        all = T)

habi.afi.funtr.iucn.syn <-
  merge(habi.afi.funtr.iucn,
        syn.wide,
        by.x = "NomS",
        by.y  = "species",
        all = T)


#------------------------------------------#
#### match all traits in a single table ####
#------------------------------------------#

nat <- habi.afi.funtr.iucn.syn
at <- acoutr
pc <- as.data.frame(pc.coord[,1:3])
colnames(pc) <- c("label","nonphy.acou.PC1","nonphy.acou.PC2")
ppc <- phy.pca$li[,1:2]
colnames(ppc) <- c("phy.acou.PC1","phy.acou.PC2")

xt <- merge(nat,at,by.x = "NomS", by.y = 0, all = F)
xpt <- merge(xt,pc,by.x = "NomS", by.y = "label", all = F)
xppt <- merge(xpt,ppc,by.x = "NomS", by.y = 0, all = F)

write.csv2(xppt,"outputs/full_trait_table.csv")

# field descriptions for the output

meta <- rbind(
c("NomS", "scientific names"),
c("SGIo" , "Godet's Species Generalism Index"),
c("RANGE_BIRDLIFE" , "Birdlife's range size from Godet et al"),
c("Affinity index" , "Blackburn's affinity index derived from Xeno Canto"),
c("Hand-Wing.Index","from AVONET"),
c("Mass","from AVONET"),
c("Habitat","from AVONET"),
c("Trophic.Niche","from AVONET"),
c("x2019_uicn_redlist" , "IUCN's redlist status (2019)"),
c("syn1 to syn4" , "synanthropy index with CARTNAT 1st layer to 4th layer"),
c("next columns", "acoustic traits (refer to dedicated trait table)"),
c("end of table" , "non phylogenetic and phylogenetic principal component axes on acoustic traits")
)
colnames(meta) <- c("column","description")

# output 

l <- list("all_traits" = xppt, "fields" = meta)
write.xlsx(l, file = "outputs/full_trait_table.xlsx")

