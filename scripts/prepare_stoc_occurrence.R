library(magrittr)


stoc <- read.csv(here::here("data", "stoc", "stoc_for_rql.csv"))
head(stoc)

col_s <- c("carre", "species", "avg.count", "occurrence", "X", "Y")

xy <- stoc[ , c("id_carre_num", "X", "Y")]
xy <- xy[!duplicated(xy$"id_carre_num"), ]

species <- sort(unique(stoc$"birdlife_sci_name"))

carres  <- sort(unique(stoc$"id_carre_num"))

dat <- expand.grid(carres, species)
colnames(dat) <- c("id_carre_num", "birdlife_sci_name") 


dat <- merge(dat, stoc[ , c("id_carre_num", "birdlife_sci_name", "mean_abundance")], 
             by = c("id_carre_num", "birdlife_sci_name"), all = TRUE)

dat$"mean_abundance" <- ifelse(is.na(dat$"mean_abundance"), 0, dat$"mean_abundance")


dat <- merge(xy, dat, by = "id_carre_num", all = TRUE)

dat$"occurrence" <- 0
dat[which(dat$"mean_abundance" > 0), "occurrence"] <- 1

dat <- dat[ , c(1, 4:6, 2:3)]

colnames(dat) <- col_s

write.csv(dat, here::here("outputs", "stoc_occurrences_at_square_level.csv"), row.names = FALSE)

# stoc.yr <- read.csv(here::here("data", "stoc", "stoc_yearly_abundances.csv"))
# 
# stoc.xy <- stoc[ , c("id_carre", "X", "Y")]
# stoc.xy <- stoc.xy[!duplicated(stoc.xy$"id_carre"), ]
# 
# 
# # add missing absences
# rg.yr        <- unique(stoc.yr$"year")
# carre.list   <- unique(stoc.yr$"id_carre")
# species.list <- unique(stoc.yr$"birdlife_sci_name")
# 
# stoc.index <- stoc.yr[ , c("id_carre","year","birdlife_sci_name")]
# 
# stoc.index.all <- stoc.index %>% 
#   dplyr::count(id_carre, year, birdlife_sci_name) %>% 
#   tidyr::complete(id_carre, year, birdlife_sci_name, fill = list(n = 0))
# 
# stoc.all  <- merge(stoc.yr, stoc.index.all, 
#                    by = c("id_carre", "year", "birdlife_sci_name"), all = TRUE)
# 
# stoc.all2 <- stoc.all[ , c("id_carre", "year", "birdlife_sci_name", "abundance")]
# stoc.all2[is.na(stoc.all2$"abundance"), "abundance"] <- 0
# 
# 
# # average abundance over years 
# stoc.avg <- aggregate(stoc.all2$"abundance",
#     by = list(stoc.all2$"id_carre", stoc.all2$"birdlife_sci_name"),
#     FUN = "mean"
#   )
# 
# colnames(stoc.avg) <- c("carre", "species", "avg.count")
# stoc.avg$"occurrence" <- 0
# stoc.avg[which(stoc.avg$"avg.count" > 0), "occurrence"] <- 1
# 
# 
# # add coordinates
# stoc.avg <- merge(stoc.avg, stoc.xy, by.x = "carre", by.y = "id_carre", all = TRUE)

# write.csv(stoc.avg, here::here("outputs", "stoc_occurrences_at_square_level.csv"), row.names = FALSE)
