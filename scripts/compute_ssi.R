source(here::here("functions", "ssi.R"))

clock <- Sys.time()

stoc <- read.csv(here::here("outputs", "stoc_occurrences_at_square_level.csv"))

stoc <- stoc[stoc$"occurrence" > 0, ]

# source(here::here("analyses", "species_list.R"))

# stoc <- stoc[which(stoc$"species" %in%
#                      c(common.species, uncommon.species, rare.species)), ]

stoc_sf <- sf::st_as_sf(stoc, coords = c("X", "Y"), crs = "epsg:27572")

stoc_sf <- sf::st_transform(stoc_sf, crs = "epsg:2154")

stoc <- data.frame(sf::st_drop_geometry(stoc_sf), sf::st_coordinates(stoc_sf))
colnames(stoc)[1] = "Cell"
colnames(stoc)[2] = "Species"
colnames(stoc)[4] = "Abundance"


layers <- 1:4

for (i in layers) {
  
  layer <- terra::rast(here::here("data", "cartnat", 
                                  paste0("CARTNAT_LAYER", i, 
                                         "_250M_BILINEAR.tif")))
  
  ssi_res <- ssi(r = layer, data = stoc, resolution = 1, sim = 99,
                 threshold = 30, n_threads = 20)

  save(ssi_res, file = here::here("outputs", paste0("SSI_CARTNAT_LAYER_", i)))
}

print(Sys.time() - clock)
