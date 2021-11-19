library(foreach)
library(doParallel)
library(spatstat.data)
library(nlme)
library(rpart)
library(spatstat)

calcEnv <- function(A,fname){

  boxr    <- box3(xrange = c(1,2500),yrange = c(1,2500),zrange = c(1,917),unitname = "microns") #1,max(A[,3])
  coorA   <- pp3(A$V1,A$V2,A$V3,boxr,marks = NULL)
#  calcEnv <- envelope.pp3(coorA, K3est, correction = c("isotropic"), nsim=4, nrank=1)#, transform = expression(. -4*pi*r^3/3))
  calcEnv <- K3est(coorA, correction = c("isotropic"))
  write.csv(calcEnv,paste("../resEnv_CA_s12/Kenv_CA_series_12",fname,"csv",sep = "."),row.names = FALSE)
  return(fname)
}

# Series 1
dat.files     <- list.files(path="../res_coord_sim_series_12/run4_adhes-5_phenchange0_set12_best_params/scale_coord/",recursive=T,pattern="*.txt",full.names=T) # Must be scaled coordinates!!!
coord         <- lapply(dat.files, read.delim, sep=",", header=FALSE)
names(coord)  <- dat.files

# Use the environment variable SLURM_CPUS_PER_TASK to set the number of cores.
# This is for SLURM. Replace SLURM_CPUS_PER_TASK by the proper variable for your system.
# Avoid manually setting a number of cores.
ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores) # Shows the number of Parallel Workers to be used
print(ncores)                    # this how many cores are available, and how many you have requested.
getDoParWorkers()                # you can compare with the number of actual workers

nsamples <- 84

# be careful! foreach() and %dopar% must be on the same line!
foreach(i=1:nsamples, .verbose=T) %dopar% {
  A<-dat.files[i]
  fname<-basename(A)
  fname
  dat<-coord[[A]]
  calcEnv(dat,fname)
}

#################### END #################################################
##########################################################################

