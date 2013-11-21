####################################
# Author: Jeremy Werner            #
# Assignment: STA 205 HW2          #
# Due date: November 13, 2013      #
####################################

rm(list = ls()) #Starting with a clean workspace (for safety)

mini = FALSE

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

###################
sim_start <- 1000 #lowest possible dataset number
s_num = 5 #number of subsamples
r_num = 50 #number of boostraps per subsample
gamma = .7
###################

if(length(args)==0){
	job_num = 1
    sim_num = sim_start + job_num
    set.seed(121231)
}else{
    job_num = as.numeric(args[1])
    sim_num <- sim_start + job_num
	# Find r and s indices: 
	s = ((job_num-1) %% s_num) + 1
	r = ceiling(job_num/s_num)
	sim_seed = (762*s + 121231)
	
	#Setting seed to get the correct rows from the sample function based on the value of s
	set.seed(sim_seed)
}

cat(paste("\nAnalyzing simulation number ",sim_num,"...\n\n",sep=""))
cat(paste("S index:", s, "   R index:", r, "\n\n"))

#============================== Run the simulation study ==============================#

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

#datapath = "/Users/jmwerner1123/Dropbox/GitHub/Invisible/STA250/HW2/BLB"
#outpath = "/Users/jmwerner1123/Dropbox/GitHub/Invisible/STA250/HW2/BLB/output"

# mini or full?
if(mini){
	rootfilename <- "blb_lin_reg_mini"
}else{
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
data_file = paste0(datapath,"/",rootfilename, ".desc")

# Attach big.matrix :
data = attach.big.matrix(data_file)

# Remaining BLB specs:
dims = dim(data)
n = dims[1] #Number of observations
d = dims[2] - 1 #Number of covariates
b = ceiling(n ^ gamma) #Size of bootstrap dataset

# Extract the subset:
subset_index = sample(1:n, size = b, replace = FALSE)
subset = data[subset_index,]

y = subset[,dims[2]] #Response
sub_data = subset[,-dims[2]]

# Reset simulation seed:
set.seed(762*sim_num + 121231) #unique seed for each job

# Bootstrap dataset:
bootstrap_data_weights = rmultinom(1, n, rep(1/b,b))

# Fit lm:
linear_model = lm(y~(0+sub_data), weights = bootstrap_data_weights)

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s),"_",sprintf("%02d",r),".txt")

# Save estimates to file:
write.table(linear_model$coeff, file = outfile, sep = ",", 
	col.names = TRUE, row.names = FALSE)


cat("Program complete, self destruct sequence initiated... \n")
q("no")
