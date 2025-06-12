
library("EnvRtype")

writeLines(" - GET DATA")
data("maizeYield") # toy set of phenotype data (grain yield per environment)
data("maizeG"    ) # toy set of genomic relationship for additive effects 
data("maizeWTH")   # toy set of environmental data
y   = "value"      # name of the vector of phenotypes
gid = "gid"        # name of the vector of genotypes
env = "env"        # name of the vector of environments
maizeYield <- droplevels(maizeYield) # make sure the nlevels of GID matchs with the G matrix dimension

writeLines(" - GET KERNELS")
W  = W_matrix(env.data = maizeWTH, var.id = c("FRUE",'PETP',"SRAD","T2M_MAX"), statistic = 'mean')
## KG and KE might be a list of kernels
KE = list(W = env_kernel(env.data = W)[[2]])
KG = list(G=maizeG);

writeLines(" - CREATE KERNEL MODELS")
## Creating kernel models with get_kernel
MM    = get_kernel(K_G = KG, y = y, gid = gid, env = env, data = maizeYield,model = "MM") 
MDs   = get_kernel(K_G = KG, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs") 
EMM   = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMM") 
EMDs  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMDs") 
RMMM  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMM") 
RNMDs = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 


writeLines(" - FIT KERNEL MODELS")
fixed = model.matrix(~0+env, maizeYield)
MDs   = get_kernel(K_G = KG, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs") 
fit   = kernel_model(y = y,env = env,gid = gid, data = maizeYield,random = MDs,fixed = fixed)

print(fit$yHat)    # predicted phenotype values

print("Variance Components")
print(fit$VarComp) # variance components and confidence intervals

print("full model output")
print(fit$BGGE)    # full output of Hierarchical Bayesian Modeling


