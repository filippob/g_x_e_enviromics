library("readr")
library("pedigree")
library("tidyverse")
library("data.table")


# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/ciampolini/anafibj',
    pheno = 'Analysis/prepped_files/pheno.dat',
    ped = 'data/ped.csv.gz', ## pedigree file or file with genetic marker data
    outdir = 'Analysis/prepped_files',
    force_overwrite = FALSE
  ))
  
}

## 1. PHENOTYPES
fname = file.path(config$base_folder, config$pheno)
pheno = fread(fname)

## 2. PEDIGREE / GENETIC MARKERS
fname = file.path(config$base_folder, config$ped)
geno = fread(fname)
names(geno) <- c("id", "sire", "dam", "sex", "birthdate")

print(paste("N. of animals with pedigree data:", sum(pheno$sire %in% geno$id)))


## MAKE KINSHIP MATRIX
writeLines(" - making the kinship matrix")
animals = ifelse(geno$id %in% pheno$sire, TRUE, FALSE)
sum(animals)

ped <- geno[,1:3]
ped <- ped |> mutate(across(everything(), as.factor))
str(ped)

summary(countGen(ped))
vec = trimPed(ped = ped, data = animals, ngenback = 15)
sum(vec)
ped <- ped[vec,]

ped <- add.Inds(ped)
ped[is.na(ped)] <- 0

ped <- ped |>
  filter(id != 0)

library("AGHmatrix")

data(ped.mrode)
K <- Amatrix(ped.mrode, ploidy=2)

K <- Amatrix(data = ped, ploidy=2, dominance = FALSE)
vec = colnames(K) %in% pheno$sire
kk = K[vec,vec]

fname = file.path(config$base_folder, config$outdir, 'kinship.png')

png(fname, width = 850, height = 850, units = "px", res = 250)
heatmap(kk)
dev.off()

fname = file.path(config$base_folder, config$outdir, 'kinship.RData')
save(kk, file = fname)

