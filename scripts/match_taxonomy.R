#---------------------------------------------------------------#
##### match taxonomy in the 10 CD's for trait extraction #####
# this script is aimed to create a taxonomic index for the 
# names of the files in the 10 CDs and scientific / english names
# J.-Y. Barnagaud - may 24
#---------------------------------------------------------------#

library(readxl)
library(stringr)
library(dplyr)

## data

cd.files <- read_excel("data/list_files_10_cds.xlsx", sheet = 1)
lof20 <- read_excel("data/LOF2020IOC10-1012020.xlsx", sheet = 1)

lof20.na <- subset(lof20,!is.na(`Nom Français` ))

## match names

# names in cd list

name.tempo <- gsub( ".mp3","", cd.files$file) 
cd.files$name <- gsub("[[:digit:]]", "", name.tempo) 
cd.files$name <- trimws(cd.files$name,"left")
cd.files.df <- as.data.frame(cd.files)

# names in lof

lof20.na$latin.name <- word(lof20.na$`Nom scientifique / Autorité`, 1,2, sep=" ")
lof20.na.df <- as.data.frame(lof20.na)
lof20.na.df$`Nom Français` <- tolower(lof20.na.df$`Nom Français`)

## remove accents in both tables (otherwise issues with matching)

cd.files.df$name <- tolower(cd.files.df$name)
cd.files.df2 <-  cd.files.df %>%
  mutate(name2 = stringi::stri_trans_general(str = name, id = "Latin-ASCII"))


lof20.na.df2 <-  lof20.na.df %>%
  mutate(fr.name = stringi::stri_trans_general(str = `Nom Français`, id = "Latin-ASCII"))


## join tables

cd.allnames <- merge(cd.files.df2,lof20.na.df2,by.x= "name2",by.y = "fr.name",all.x=T,all.y = F)

## missing names

# some names have been manually updated in the below file because they were missing
# and not easily retrievable automatically due to various issues: 

mis.names <- read.csv2("data/missing_latin_names_10cd.csv")

cd.allnames2 <- cd.allnames
cd.allnames2[which(cd.allnames2$name2%in%mis.names$name2),"latin.name"] <- mis.names$latin.name
  
## export

cd.index <- cd.allnames2[,c("path","file","name","Nom scientifique / Autorité","latin.name")]
write.csv(cd.index,file = "data/index_10cd_latin.csv")
