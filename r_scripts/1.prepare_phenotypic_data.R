library("readr")
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
    y = 'data/feno_pisa.dat',
    env = 'data/herd_codici_pisa.csv',
    ped = 'data/ped.csv.gz',
    outdir = 'Analysis/prepped_files',
    force_overwrite = FALSE
  ))
  
}


## 1. PHENOTYPES

fname = file.path(config$base_folder, config$y)

## read_fwf: this is more efficient, but needs some tweaking with the code
# feno = read_fwf(fname, fwf_widths(c(8, 5, 7, 4, 2, 2, 3, 3, 1, 1, 1, 1, 3, 4, 4, 5, 2, 2, 3, 1, 4, 5, 5, 5), 
#                            c("id", "skip", "herd", "rec_year", "rec_mo", "rec_day", "age", "DIM", "nm", "ese", "ncon", 
#                              "ncam", "milk", "fat", "prot", "scc", "lactose", "skip2","dp", "cl", "flag", 
#                              "THI", "humidity", "max_temp")))

## this comes from Denis Laloe 1.lire.Rmd
feno = read.fortran(fname,format=c("A8","5X","A7","F4","2F2","F3","F3","4A1","F3","F4","F4","A5","A2","2X","A3","A1","A4","A5","A5","A5"))
names(feno) <- c("Id","Herd","rec_year","rec_month","rec_day","Age","DIM","Nm","Ese","Ncon","Ncam","Milk","Fat","Prot","scc","lactose","Dp","Classe","Flags","THI","Hum","MaxT")

feno <- feno |>
  mutate(scc = as.numeric(str_trim(scc, "both")),
         Herd = str_trim(Herd, "both"),
         THI = as.numeric(str_trim(THI, "both")),
         Hum = as.numeric(str_trim(Hum, "both")),
         MaxT = as.numeric(str_trim(MaxT, "both")),
         Dp = as.numeric(str_trim(Dp, "both"))
  )

summary(feno$scc)

## 2. PEDIGREE

fname = file.path(config$base_folder, config$ped)
ped = fread(fname)
names(ped) <- c("id", "sire", "dam", "sex", "birthdate")

print(paste("N. of phenotypic records is:", nrow(feno)))
ids = unique(feno$Id)

print(paste("N. of unique cow IDs is:", length(ids)))

print(paste("This means that the averege n. of records per cow is:", round(nrow(feno)/length(ids)),2))


if (sum(ids %in% ped$id) == length(ids)) {
  
  print("ALL cows with phenoypic records have also a pedigree record")
} else print("NOT all cows with phenoypic records have also a pedigree record")


## 3. HERDS (ENVIRONMENTS)

fname = file.path(config$base_folder, config$env)

env = fread(fname, header = TRUE, fill = TRUE)

## removing records with missing data (hack to make it work: better way to read file with unequal record length needed: 
## here we have a couple of records with a missing value represented by one additional space)
env[is.na(env$rad_long),]
env <- na.omit(env)


herds = unique(feno$Herd)
sum(herds %in% env$herd) == length(herds) ## two herds are missing, since we removed records with NA

herds <- select(env, c(herd,prov)) |>
  unique()

feno$prov = herds$prov[match(feno$Herd, herds$herd)]
feno$sire = ped$sire[match(feno$Id, ped$id)]

#####################
### select and filter
#####################

feno_red <- feno |>
  select(Id, prov, rec_year, Milk, DIM, Age, THI, Hum, MaxT, sire, Dp, rec_month, rec_day)

feno_red <- na.omit(feno_red)
feno_red <- filter(feno_red, sire > 0)

sire_rec <- feno_red |>
  group_by(sire, prov, rec_year) |>
  summarise(N = n(), MY = mean(Milk))

temp <- sire_rec |>
  group_by(prov, rec_year) |>
  summarise(tot = sum(N)) |>
  spread(key = rec_year, value = tot)

fwrite(x = temp, file = "aspa_2025/records.csv")

## from summary stats: we select CN (Piemonte), RM (Lazio), MN (est Lombardia); year 2021
temp <- feno_red |>
  select(sire,prov,rec_year) |>
  filter(rec_year == 2021) |>
  select(-rec_year) |>
  group_by(sire, prov) |>
  summarise(N= n()) |>
  spread(key = prov, value = N) |>
  arrange(desc(CR))

M <- temp[,-1]
n_recs = rowSums(M, na.rm = TRUE)
M <- !is.na(M) 
n_provs = rowSums(M)
df <- data.frame(sire = temp$sire, n_rec = n_recs, n_prov = n_provs)
df <- df |> arrange(desc(n_prov))

tot_provs = ncol(temp)-1
df$ratio = df$n_prov/tot_provs

temp <- feno_red |>
  select(sire,prov,rec_year) |>
  filter(rec_year == 2021) |>
  select(-rec_year) |>
  group_by(sire, prov) |>
  summarise(N= n()) |>
  spread(key = sire, value = N)

M <- temp[,-1]
n_recs = rowSums(M, na.rm = TRUE)
M <- !is.na(M) 
n_sires = rowSums(M)
df <- data.frame(prov = temp$prov, n_rec = n_recs, n_sire = n_sires)
df <- df |> arrange(desc(n_sire))

tot_sires = ncol(temp)-1
df$ratio = df$n_sire/tot_sires

temp <- feno_red |>
  filter(prov %in% c("CR","MN","CN"), rec_year == 2021)

X <- temp |>
  group_by(sire, prov) |>
  summarise(DIM = mean(DIM), age = mean(Age), dp = mean(Dp), 
            rec_month = mean(rec_month), rec_day = mean(rec_day))

Y <- temp |>
  select(-c(Id)) |>
  group_by(sire, prov) |>
  summarise(milk = mean(Milk), age = mean(Age), DIM = mean(DIM), 
            THI = mean(THI), maxT = mean(MaxT), humidity = mean(Hum))

Y <- Y |>
  select(-c("DIM","age","THI","humidity","maxT")) |>
  spread(key = prov, value = milk)

## creating output folder if not there
dir.create(file.path(config$base_folder, config$outdir), showWarnings = FALSE, recursive = TRUE)

fname = file.path(config$base_folder, config$outdir, 'pheno.csv')
fwrite(x = Y, file = fname, col.names = TRUE)

X <- Y |>
  gather(key = "prov", value = "milk", -sire) |>
  left_join(X, by = c("sire", "prov")) |>
  select(-milk)

fname = file.path(config$base_folder, config$outdir, 'covariates.csv')
fwrite(x = X, file = fname, col.names = TRUE)
