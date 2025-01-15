######################################################################################################################

# GxE using marker-by-environment interactions


# (2) Using genomic relationships

#A model equivalent to the one presented above can be implemented using G-matrices (or factorizations of it) 
# with off-diagonal blocks zeroed out for interactions

######################################################################################################################

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
    base_folder = '~/Documents/ciampolini/anafibj/',
    y = 'Analysis/prepped_files/pheno.csv',
    X = 'Analysis/prepped_files/covariates.csv',
    K = 'Analysis/prepped_files/kinship.RData',
    test_prop = 0.2,
    evd_threshold = 1e-3,
    nIter = 10000,
    burnIn = 2000,
    thin = 5, ## default value in BGLR is 5
    outdir = 'Analysis/BGLR',
    prefix = "GxE_kin_",
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
library("BGLR")
library("knitr")
library("BGData")
library("tidyverse")
library("data.table")

# READ THE DATA -------------------------------------------------------------------
## phenotypes (same phenotypes, multiple envs (maybe also multiple years))
fname = file.path(config$base_folder, config$y) 
Y = fread(fname, header = TRUE) # grain yield evaluated in 4 different environments (599 samples, 4 columns)

## covariates
fname = file.path(config$base_folder, config$X) 
X <- fread(fname) ## 599 samples, 1279 markers (DArT markers --> 0/1)
dim(X)

## kinship
fname = file.path(config$base_folder, config$K) 
load(fname) ## 599 samples, 1279 markers (DArT markers --> 0/1)
dim(kk)

ntraits = ncol(Y)-1 

print("raw correlations between phenotypes")
print(kable(round(cor(Y, use = "complete.obs"),2))) # correlation between yields in the 4 environments (round with 2 decimals)

####################
writeLines(" - preparing long vector of phenotypes")
y <- Y |> gather(key = "trait", value = "value", -sire)
y <- filter(y, !is.na(value))
v <- group_by(y, trait) |> summarise(N=n()) |> pull(N)

print(kable(v))

env = NULL
for (i in 1:length(v)) {
  print(i)
  env = c(env, rep(i,v[i]))
}

#2# Creating a Testing set
n = nrow(y)
test_size = n*(config$test_prop)
yNA <- y
tst <- sample(1:n, size=test_size, replace=FALSE)
yNA[tst, "value"] <- NA

writeLines(" - extracting eigenvalues and eigenvectors from the kinship matrix")
# EVD = eigen(kk)
# PC = EVD$vectors[,EVD$values > config$evd_threshold]   # 599 x 598 ? SELECT THRESHOLD TO RETAIN PC's
# for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }  #? eigenvectors multiplied by the sqrt of the corresponding eigenvalues

load("Analysis/BGLR/princ_comps.RData")

writeLines(" - preparing design matrices for multiple envs")
## make "p" copies of PC <-- "p" phenotypes (same pheno in "p" envs, same sample IDs: this is why it is duplicated)
row.names(PC) <- as.character(unique(X$sire))
X0 = PC[as.character(y$sire),] # Matrix for main effects 

# now interactions
for (i in 1:ntraits) {
  n = ntraits
  print(paste(i,n))
  assign(paste("X", i, sep = ""), X0)
}

#############################################
## !! MANUAL EDITING HERE IF NTRAITS > 4!! ##
#############################################
if (ntraits == 2) {
  
  for(i in 1:nrow(X0)){
    X1[i,]<-(env[i]==1)*X0[i,]
    X2[i,]<-(env[i]==2)*X0[i,]	
  }
  
  LP = list(main=list(X=X0,model='BRR'), 
            int1=list(X=X1,model='BRR'),
            int2=list(X=X2,model='BRR')
  )
} else if (ntraits == 3) {
  
  for(i in 1:nrow(X0)){
    X1[i,]<-(env[i]==1)*X0[i,]
    X2[i,]<-(env[i]==2)*X0[i,]	
    X3[i,]<-(env[i]==3)*X0[i,]
  }
  
  LP = list(main=list(X=X0,model='BRR'), 
            int1=list(X=X1,model='BRR'),
            int2=list(X=X2,model='BRR'),
            int3=list(X=X3,model='BRR')
  )
} else if (ntraits == 4) {
  
  for(i in 1:nrow(X0)){
    X1[i,]<-(env[i]==1)*X0[i,]
    X2[i,]<-(env[i]==2)*X0[i,]	
    X3[i,]<-(env[i]==3)*X0[i,]
    X4[i,]<-(env[i]==4)*X0[i,]
  }
  
  LP = list(main=list(X=X0,model='BRR'), 
            int1=list(X=X1,model='BRR'),
            int2=list(X=X2,model='BRR'),
            int3=list(X=X3,model='BRR'),
            int4=list(X=X4,model='BRR')
  )
}

###########################
# BGLR MODEL -------------------------------------------------------------------
print("Running the BGLR model - kinship matrix")
experiment = "milk"
dirname = file.path(config$base_folder,config$outdir,experiment)
dir.create(dirname, recursive = TRUE, showWarnings = FALSE)
outpath = file.path(dirname, config$prefix)
fmGRM = BGLR(y=yNA$value,ETA=LP,
             nIter=config$nIter, burnIn=config$burnIn, thin = config$thin,
             saveAt=outpath, groups=env)

print("Writing out results")
fname = paste(config$prefix, "BRR_res.RData", sep="")
to_save = list("fmGRM"=fmGRM, "y"=y, "test_rows" = tst)
save(to_save, file = file.path(dirname, fname))

print("DONE!")



