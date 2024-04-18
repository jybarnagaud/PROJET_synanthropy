#---------------------------------------------------------------------------#
#### coinertia analysis between acoustic and non-acoustic traits #####
          # Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
                        # CESAB/ACOUCENE project
            # subproject on acoustic vs non acoustic traits
                              # april 2024
#---------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(viridis)
library(ade4)
library(adegraphics)
library(factoextra)
library(PCAtest)
library(PerformanceAnalytics)

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
  traits[, c(
    "duration_song",
    "centroid_f",
    "LFC_spectral_cover",
    "MFC_spectral_cover",
    "HFC_spectral_cover",
    "Hf_spectral_entropy",
    "number_of_freq_peaks",
    "syllable_duration_mean",
    "sound_per_vocalization_ratio"
  )]

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


# gap filling (mean) for missing values

geo.traits.filled <- geo.traits %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE))))

fun.traits.filled <- fun.traits %>%
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

## single - table ordinations --------------------------------------------------
