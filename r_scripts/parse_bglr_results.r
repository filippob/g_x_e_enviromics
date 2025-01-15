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
    # base_folder = '~/Documents/freeclimb/g_x_e',
    base_folder = '~/Documents/ciampolini/anafibj/',
    X = 'Analysis/prepped_files/covariates.csv',
    result_folder = 'Analysis/BGLR',
    outdir = 'Analysis/BGLR/results',
    traits = 'milk',
    prefix = "GxE_kin_",
    ntraits = 3,
    burnIn = 100,
    thin = 5, ## default thin in BGLR is 5
    force_overwrite = FALSE
  ))
  
}

library("stringr")
library("data.table")

traits = strsplit(config$traits, ",")[[1]]

res = data.frame("trait"=NULL, "test_corr"=NULL, "train_corr"=NULL, 
                 "test_rmse" = NULL, "train_rmse" = NULL,
                 "h2_1"=NULL, "h2_2"=NULL, "h2_3"=NULL,
                 "rg1_2"=NULL, "rg1_3"=NULL)
for (trt in traits) {
  
  trt = str_trim(trt, side = "both")
  fpath = file.path(config$base_folder, config$result_folder, trt)

  filenames = list.files(fpath, pattern="*.RData", full.names=TRUE)
  fname = filenames[grepl(pattern = config$prefix, filenames)][1]
  print(paste("Loading data from",fname))
  load(fname)
  
  y <- to_save$y
  fmGRM <- to_save$fmGRM
  tst = to_save$test_rows
  
  #5# Assesment of correlation in TRN and TST data sets
  test_corr = cor(fmGRM$yHat[tst],y$value[tst])
  train_corr = cor(fmGRM$yHat[-tst],y$value[-tst]) #TRN
  
  test_rmse = sqrt(sum((fmGRM$yHat[tst]-y$value[tst])^2)/length(tst))
  train_rmse = sqrt(sum((fmGRM$yHat[-tst]-y$value[-tst])^2)/length(y$value[-tst]))
  
  ## residual variance
  filenames = list.files(fpath, pattern="*varE.dat", full.names=TRUE)
  fname = filenames[grepl(pattern = config$prefix, filenames)][1]
  print(paste("varE:",fname))
  varE = fread(fname)
  varE = as.data.frame(varE[-c(1:(config$burnIn/config$thin)), ])
  
  ## main genetic variance
  filenames = list.files(fpath, pattern="*main_varB.dat", full.names=TRUE)
  fname = filenames[grepl(pattern = config$prefix, filenames)][1]
  print(paste("varU main:",fname))
  varU_main = fread(fname)
  varU_main = varU_main$V1[-c(1:(config$burnIn/config$thin))]
  
  ## GxE variances
  filenames = list.files(fpath, pattern="*int[0-9]_varB.dat", full.names=TRUE)
  filenames = filenames[grepl(pattern = config$prefix, filenames)]
  print(filenames)

  for (i in 1:length(filenames)) {
    
    fname = paste(config$prefix, 'ETA_int', i, '_varB.dat', sep="")
    temp = scan(file.path(fpath, fname))[-c(1:(config$burnIn/config$thin))]
    assign(paste("varU_int", i, sep = ""), temp)
    assign(paste("varU", i, sep=""), get(paste("varU_int", i, sep = "")) + varU_main)
    assign(paste("h2", i, sep="_"), get(paste("varU", i, sep=""))/(get(paste("varU", i, sep=""))+varE[,i]))
    
    
  }
  
  ## pairwise combinations
  m = combn(length(filenames),2)
  lcomb = split(m, rep(1:ncol(m), each = nrow(m)))
  
  ## GENETIC CORRELATIONS
  for (i in 1:length(lcomb)) {
    
    tmp = varU_main/sqrt(get(paste("varU",lcomb[[i]][1],sep="")) * get(paste("varU",lcomb[[i]][2],sep="")))
    assign(paste("rg",lcomb[[i]][1],lcomb[[i]][2],sep="_"), tmp)
  }
  
  temp = data.frame("trait"=trt, "test_corr"=test_corr, "train_corr"=train_corr,
                    "test_rmse"=test_rmse, "train_rmse"=train_rmse,
                    "h2_1" = ifelse(exists("h2_1"),mean(h2_1),NA),
                    "h2_2" = ifelse(exists("h2_2"),mean(h2_2),NA),
                    "h2_3" = ifelse(exists("h2_3"),mean(h2_3),NA),
                    "rg1_2" = ifelse(exists("rg_1_2"),mean(rg_1_2),NA),
                    "rg1_3" = ifelse(exists("rg_1_3"),mean(rg_1_3),NA),
                    "rg2_3" = ifelse(exists("rg_2_3"),mean(rg_2_3),NA)
                    )
  res = rbind.data.frame(res,temp)
  
  rm(y,fmGRM,tst,m)
  rm(list=ls(pattern="var"))
  rm(list=ls(pattern="h2"))
  rm(list=ls(pattern="rg"))
}

dir.create(file.path(config$base_folder, config$outdir), showWarnings = FALSE)
fname = file.path(config$base_folder, config$outdir, "results.csv")
fwrite(x = res, file = fname, sep=",")

print("DONE!!")



