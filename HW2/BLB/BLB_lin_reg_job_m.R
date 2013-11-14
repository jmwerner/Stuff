# Fit simple linear regression on d=1,000 covariates with *no* intercept
# Errors ~ N(0,sigma^2)
# Goal is to find SE(beta.hat.1),...,SE(beta.hat.1000)

# Data description:
# blb_lin_reg_data.txt: n=1,000,000 observations. Each row corresponds to an observation,
# each column to a variable. First d=1,000 columns are covariates (X1,...,X1000), 
# final column corresponds to the response variable y
# Mini dataset d=40 covariates and n=10,000 observations

# Clear out everything in memory for a nice clean run and easier debugging
rm(list=ls())

s = 5           # s = number of distinct subsamples 
r = 50          # r = number of bootstrap samples to do for each distinct subsample
mini <- FALSE   # use the mini dataset or the full dataset?
# args = 18
#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])

  # Get the job number from the argument 1:(s*r) (i.e. 1:5*50 = 1:250)
  job.num = as.numeric(args[1])
  
  # Get the s_index by using mod s.  
  # Also, if job.num mod s == 0, then it is subsample dataset s (i.e. 5)
  s_index = job.num %% s
  if (s_index == 0){
    s_index = s
  }
  
  # Get the r_index by using ceiling(job.num/s).  
  # Also, if job.num mod r == 0, then it is bootstrap sample r (i.e. 50) within subsample dataset s
  r_index = ceiling(job.num / s)
  if (r_index == 0){
    r_index = r
  }
  
  # The first seed must be a function of the s_index to ensure that the subsample is the same
  # for same values of the s_index
  sim_seed <- (762*(s_index) + 121231)
  set.seed(sim_seed)
}

# Some checks that go into the .out file
cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))
cat(paste("\nRunning s_index ",s_index," r_index ",r_index," seed ",sim_seed," job.num ",job.num,"...\n\n",sep=""))


#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
# datapath = "C:/Users/Michael/Documents/Michael UC Davis/STA 250 Adv Stat Computing/HW2/"
datapath <- "/home/pdbaines/data/"
outpath <- "output/"


# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
datafile = paste0(datapath, rootfilename, ".desc")

# Set up I/O stuff:

# Attach big.matrix :
data.full = attach.big.matrix(datafile)

# Remaining BLB specs:
n = nrow(data.full)     # n = nrows of the full dataset
d = dim(data.full)[2]-1 # d = number of covariates
gamma = 0.7
b = ceiling(n^gamma)    # b = size of each subset b<n, taken without replacement from the full dataset
                        # then resample n points from b<n

# Extract the subset:
# Get b row indices from 1:n without replacement to use for each of the s subsamples
# based on the seed set above and use the same seed for same values of s_index
row.indices = sample(1:n, size=b, replace=FALSE)
X.sub.sample = data.full[row.indices,1:d]
Y.sub.sample = data.full[row.indices,d+1]

# Checks the subsets:
# dim(X.sub.sample) 
# dim(Y.sub.sample) 
outfile = paste0("dump/","BLB_sample_indices_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")
write.table(x=cbind(sim_seed,row.indices),file=outfile, sep=",", col.names=TRUE, row.names=FALSE)

# Reset simulation seed:
# The seed for the bootstrap sample must be different for each bootstrap sample r within
# each subsample s, so make the seed a function of s_index and r_index with no overlap
sim_seed <- (762*(s_index) + 121231 + r_index)
set.seed(sim_seed)

# Bootstrap dataset:
# Select bootstrap sample of size n with replacement from subsample
# Do this using rmultinom to distribute n balls into b boxes
# where each box has equal probability, i.e. 1/b
# This gives us the weights (i.e. # of replications) to put on each of the b row indices in 
# our subsample so that we get a bootstrap sample of size n
data.weights = rmultinom(1, size = n, rep(1/b, b))

# Check the weights
outfile = paste0("dump/","BLB_bootstrap_weights_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")
write.table(x=data.weights,file=outfile, sep=",", col.names=TRUE, row.names=FALSE)

# Fit the linear regression using lm and get the coefficients:
model = lm(Y.sub.sample ~ 0 + X.sub.sample, weights=data.weights)
beta.hat = model$coefficients

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:
write.table(x=beta.hat,file=outfile, sep=",", col.names=TRUE, row.names=FALSE)

cat("done. :)\n")

q("no")
