###########################
# Author: Jeremy Werner	  #
# STA 250 Driver for      #
# Monte Carlo Integration #
# with CUDA               #
###########################

require(RCUDA)

##########################
## Function definitions ##
##########################

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

"stnorm_pdf" = function(x){
	return(exp(-(x^2)/2)/sqrt(2*pi))	
}# stnorm_pdf

"normal_table_mc_cpu" = function(max_quantile, highest_x, N){
	table_quantiles = (1:(max_quantile * 100)) / 100

	n = length(table_quantiles)
	normal_vector = rep(0, n)
	
	lo_y = 0
	hi_x = highest_x
	hi_y = 1
	
	for(iter in 1:n){
		lo_x = table_quantiles[iter]
		x = runif(N, min = lo_x, max = hi_x)
		y = runif(N, min = lo_y, max = hi_y)
		vect = y < stnorm_pdf(x)
		normal_vector[iter] = sum(vect) * (hi_x - lo_x) * (hi_y - lo_y) / N
	}
	return(matrix(normal_vector, ncol = 10, byrow = TRUE))
}# normal_table_mc_cpu 

"normal_table_mc_vegas_cpu" = function(max_quantile, lambda, N){
	table_quantiles = (1:(max_quantile * 100)) / 100

	n = length(table_quantiles)
	normal_vector = rep(0, n)
	
	for(iter in 1:n){
		lo_x = table_quantiles[iter]
		x_new = rexp(N, lambda)
		#y = runif(N, min = lo_y, max = hi_y)
		vect = stnorm_pdf(x_new)/(lambda * exp(-lambda * x_new))
		vect[x_new < lo_x] = 0
		normal_vector[iter] = sum(vect)/ N
	}
	return(matrix(normal_vector, ncol = 10, byrow = TRUE))
}# normal_table_mc_vegas_cpu 

"normal_table_mc_gpu" = function(kernel, max_quantile, highest_x, N, verbose = FALSE){
	bg = compute_grid(N)
	grid_dims = bg$grid_dims
	block_dims = bg$block_dims
	
	nthreads = prod(grid_dims) * prod(block_dims)
	if(verbose){
		cat("Total number of threads to be launched = ",nthreads,"\n")
	}
	if(nthreads < N){
		stop("Grid is not large enough! Self destruct sequence initiated...")
	}
	
	table_quantiles = (1:(max_quantile * 100)) / 100

	n = length(table_quantiles)
	normal_vector = rep(0, n)
	
	lo_y = 0
	hi_x = highest_x
	hi_y = 1
	vect = rep(0L,N)
	
	for(iter in 1:n){
		lo_x = table_quantiles[iter]

		vals_dev = copyToDevice(vect)
		
		.cuda(kernel, vals_dev, N, N, lo_x, hi_x, gridDim = grid_dims, blockDim = block_dims)
		
		vect = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="int")
		
		normal_vector[iter] = sum(vect) * (hi_x - lo_x) * (hi_y - lo_y) / N
	}
	return(matrix(normal_vector, ncol = 10, byrow = TRUE))
}# normal_table_mc_gpu

"normal_table_mc_vegas_gpu" = function(kernel, max_quantile, N , lambda, verbose = FALSE){
	bg = compute_grid(N)
	grid_dims = bg$grid_dims
	block_dims = bg$block_dims
	
	nthreads = prod(grid_dims) * prod(block_dims)
	if(verbose){
		cat("Total number of threads to be launched = ",nthreads,"\n")
	}
	if(nthreads < N){
		stop("Grid is not large enough! Self destruct sequence initiated...")
	}
	
	table_quantiles = (1:(max_quantile * 100)) / 100

	n = length(table_quantiles)
	normal_vector = rep(0, n)

	vect = rep(0.0,N)
	
	for(iter in 1:n){
		quan = table_quantiles[iter]

		vals_dev = copyToDevice(vect)
		
		.cuda(kernel, vals_dev, N, N, quan, lambda, gridDim = grid_dims, blockDim = block_dims)
		
		vect = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="int")
		
		normal_vector[iter] = sum(vect) / N
	}
	return(matrix(normal_vector, ncol = 10, byrow = TRUE))
}# normal_table_mc_gpu


"normal_table_mc_gpu_2" = function(kernel, max_quantile, highest_x, N, verbose = FALSE){
	table_quantiles = (1:(max_quantile * 100)) / 100
	n = length(table_quantiles)
	
	bg = compute_grid(n)
	grid_dims = bg$grid_dims
	block_dims = bg$block_dims
	
	nthreads = prod(grid_dims) * prod(block_dims)
	if(verbose){
		cat("Total number of threads to be launched = ",nthreads,"\n")
	}
	if(nthreads < n){
		stop("Grid is not large enough! Self destruct sequence initiated...")
	}

	lo_y = 0
	hi_x = highest_x
	hi_y = 1
	vect = rep(0.0,n)
	
	vals_dev = copyToDevice(vect)
		
	.cuda(kernel, vals_dev, n, N, table_quantiles, n, highest_x, gridDim = grid_dims, blockDim = block_dims)
		
	vect = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="int")


	return(matrix(vect, ncol = 10, byrow = TRUE))
}# normal_table_mc_gpu_2

########################################################################################

################
## Main Block ##
################

# GPU SETUP (truncnorm kernel developed in question 1)
cat("Setting cuGetContext(TRUE)...\n")
cuGetContext(TRUE)
cat("done. Profiling CUDA code...\n")

cat("Loading module...\n")
m = loadModule("mc_integration.ptx")  
cat("done. Extracting kernel...\n")
k = m$mc_integration_normal_kernel
k2 = m$mc_integration_normal_vegas_kernel
k3 = m$mc_integration_normal_kernel_2

# Inputs to be changed

k_in = 7

N_in = 10^k_in
max_x_in = 10
lambda_in = 1.5

# Execution section

####
cpu_time = system.time({
	cpu_table = normal_table_mc_cpu(max_quantile = 3, max_x_in, N = N_in)
})

gpu_time = system.time({
	gpu_table = normal_table_mc_gpu(k,max_quantile = 3, max_x_in, N = N_in)
})

cpu_vegas_time = system.time({
	cpu_vegas_table = normal_table_mc_vegas_cpu(max_quantile = 3, lambda = lambda_in, N = N_in)
})

gpu_vegas_time = system.time({
	gpu_vegas_table = normal_table_mc_vegas_gpu(k2, max_quantile = 3, N = N_in, lambda = lambda_in)
})

#gpu_time_2 = system.time({
#	gpu_table_2 = normal_table_mc_gpu_2(k3, max_quantile = 3, highest_x = max_x_in, N = N_in)
#})
####

#Output section

write.csv(cpu_table, paste0("normal_table_cpu_",k_in,".csv"))
write.csv(gpu_table, paste0("normal_table_gpu_",k_in,".csv"))
write.csv(cpu_vegas_table, paste0("normal_table_vegas_cpu_",k_in,".csv"))
write.csv(gpu_vegas_table, paste0("normal_table_vegas_gpu_",k_in,".csv"))
#write.csv(gpu_table_2, paste0("normal_table_gpu_2".csv"))

write.csv(rbind(cpu_time, gpu_time, cpu_vegas_time, gpu_vegas_time), paste0("times_",k_in,".csv"))
#print(cpu_time)
#print(gpu_time)
#print(cpu_vegas_time)
#print(gpu_vegas_time)
#print(gpu_time_2)

q("no")

