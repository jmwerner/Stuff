####################################
# Author: Jeremy Werner            #
# Assignment: STA 205 HW4          #
# Due date: December 9, 2013       #
####################################

# Executes gibbs samping algorithm for CPU and GPU applications of a 
# Bayesian Probit regression model

# libraries
require(MASS)
require(truncnorm)
require(RCUDA)

# Command line arguments 
args <- commandArgs(trailingOnly = TRUE)
cat(paste("Command Line Arguments:", args, "\n"))

if(length(args) == 0){
	data_num = "05" #Set manually in case the user's attitude is too poor for
					#  command line arguments
}else{
	data_num = args[[1]]
}	

data = read.table(paste0("data_", data_num, ".txt"), header = TRUE)

##########################
## Function definitions ##
##########################

# Returns posterior up to proportionality
"posterior" = function(
	precision,	# (p x p) prior precision
	beta0,		# (p x 1) original beta0 vector
	X,			# (n x p) design matrix
	Z			# (n x 1) vector of Z values
){
	v = solve(precision + t(X) %*% X)
	m = v %*% (precision %*% beta0 + t(X) %*% Z)
	return(mvrnorm(1, m, v))
}# posterior

"probit_mcmc_cpu" = function(
     y,           			# vector of length n (binary)
     X,           			# (n x p) design matrix
     beta_0,      			# (p x 1) prior mean
     Sigma_0_inv, 			# (p x p) prior precision 
     niter = 2000,       		# number of post burnin iterations
     burnin = 500,	  		# number of burnin iterations
     convergence_plot = FALSE	# flag to print convergence plot
){
	N = length(y)
	p = length(beta_0)
	lowbounds = rep(0, N)
	lowbounds[!y] = -Inf
	upbounds = rep(0, N)
	upbounds[y>0] = Inf
	
	# Pre-allocate space for results
	betas_burnin = matrix(rep(0,burnin*p),nrow = burnin, ncol = p)
	betas = matrix(rep(0,niter*p),nrow = niter, ncol = p)
	
	beta_curr = as.numeric(beta_0)
	for(iter in 1:burnin){
		new_z = rtruncnorm(N, lowbounds, upbounds, mean = X %*% beta_curr, sd = rep(1,N))
		beta_curr = posterior(Sigma_0_inv, beta_0, X, new_z)
		betas_burnin[iter,] = beta_curr
	}#for	
	
	for(iter in 1:niter){
		new_z = rtruncnorm(N, lowbounds, upbounds, mean = X %*% beta_curr, sd = rep(1,N))
		beta_curr = posterior(Sigma_0_inv, beta_0, X, new_z)
		betas[iter,] = beta_curr
	}#for	
	
	if(convergence_plot){
		##Bult for 9 or less parameters, didn't feel like making the layout dynamic
		par(mfrow = c(3,3))
		lapply(1:ncol(betas), FUN = function(x, data, burn){
			plot(c(burn[i,],data[i,]), type = "l")
			abline(v = burnin, col = "RED")
		}, data = t(betas), burn = t(betas_burnin))
	}
	return(betas)
}#probit_mcmc_cpu

"compute_grid" <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
    # if...
    # N = 1,000,000
    # => 1954 blocks of 512 threads will suffice
    # => (62 x 32) grid, (512 x 1 x 1) blocks
    # Fix block dims:
    block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
    threads_per_block <- prod(block_dims)
    if (grid_nd==1){
      grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
      grid_d2 <- 1L
    } else {
      grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
      grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
    }
    grid_dims <- c(grid_d1, grid_d2, 1L)
    return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

"probit_mcmc_gpu" = function(
      y,           				# vector of length n 
      X,           				# (n x p) design matrix
      beta_0,      				# (p x 1) prior mean
      Sigma_0_inv, 				# (p x p) prior precision 
      niter,       				# number of post burnin iterations
      burnin,       			# number of burnin iterations
      verbose = FALSE,			# Print updates?
      convergence_plot = FALSE	# Plot traces of parameters?
){
	Inf_equivalent = 1e8L
	N = length(y)
	p = length(beta_0)
	lowbounds = rep(0L, N)
	lowbounds[!y] = -Inf_equivalent
	upbounds = rep(0L, N)
	upbounds[y>0] = Inf_equivalent
	
	bg = compute_grid(N)
	grid_dims = bg$grid_dims
	 block_dims = bg$block_dims
	 
	 nthreads <- prod(grid_dims)*prod(block_dims) 
	 if(verbose){
	 	cat("Total number of threads to launch = ",nthreads,"\n")
	 }
	 if (nthreads < N){
    	stop("Grid is not large enough...! Self destruct sequence initiated...")
	 }

	 
	# Pre-allocate space for results
	betas_burnin = matrix(rep(0,burnin*p),nrow = burnin, ncol = p)
	betas = matrix(rep(0,niter*p),nrow = niter, ncol = p)
	
	beta_curr = beta_0
	new_z = rep(0.0, N)
	sigma_in = rep(1,N)
	maxtries_in = 2000L
	
	vals_dev = copyToDevice(new_z)
	 sigma_dev = copyToDevice(sigma_in)
	 lo_dev = copyToDevice(lowbounds)
	 hi_dev = copyToDevice(upbounds)
		
	for(iter in 1:burnin){
		mu_in = X %*% beta_curr
		mu_dev = copyToDevice(mu_in)
		.cuda(k, vals_dev, N, mu_dev, sigma_dev, lo_dev, hi_dev, N, N, N, N, maxtries_in,
			 gridDim = grid_dims, blockDim = block_dims)
		new_z = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="float")
		beta_curr = posterior(Sigma_0_inv, beta_0, X, new_z)
		betas_burnin[iter,] = beta_curr
	}#for	
	
	for(iter in 1:niter){
		mu_in = X %*% beta_curr
		mu_dev = copyToDevice(mu_in)
		.cuda(k, vals_dev, N, mu_dev, sigma_dev, lo_dev, hi_dev, N, N, N, N, maxtries_in,
			 gridDim = grid_dims, blockDim = block_dims)
		new_z = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="float")
		beta_curr = posterior(Sigma_0_inv, beta_0, X, new_z)
		betas[iter,] = beta_curr
	}#for	
	
	if(convergence_plot){
		##Bult for 9 or less parameters, didn't feel like making the layout dynamic
		par(mfrow = c(3,3))
		lapply(1:ncol(betas), FUN = function(x, data, burn){
			plot(c(burn[i,],data[i,]), type = "l")
			abline(v = burnin, col = "RED")
		}, data = t(betas), burn = t(betas_burnin))
	}
	return(betas)
}# probit_mcmc_gpu

########################################################################################

################
## Main Block ##
################


# GPU SETUP (truncnorm kernel developed in question 1)
cat("Setting cuGetContext(TRUE)...\n")
cuGetContext(TRUE)
cat("done. Profiling CUDA code...\n")

cat("Loading module...\n")
m = loadModule("truncnorm.ptx")  
cat("done. Extracting kernel...\n")
k = m$rtruncnorm_kernel

###################################


# Define inputs
y = data[,1]
beta_0 = rep(0L,8)
X = as.matrix(data[,-1])
Sigma_0_inv = matrix(rep(0,8*8),8,8)
niter = 2000
burnin = 500

cpu_time = system.time({
	betas_cpu = probit_mcmc_cpu(y,X,beta_0,Sigma_0_inv,2000,500)
})

print(cpu_time)

gpu_time = system.time({
	betas_gpu = probit_mcmc_gpu(y, X, beta_0, Sigma_0_inv, 2000, 500)
})

print(gpu_time)

#Checking of estimates (they were bad!)
print(apply(betas_cpu,2,"mean"))
print(apply(betas_gpu,2,"mean"))

write.csv(rbind(gpu_time,cpu_time), paste0("data_", data_num,"_times.csv"))

write.csv(rbind(summary(betas_cpu), summary(betas_gpu)), paste0("data_", data_num,"_beta_summary.csv"))



###########################################################################

q("no")
