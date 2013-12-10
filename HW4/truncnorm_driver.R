require(RCUDA)
require(truncnorm)

cat("Setting cuGetContext(TRUE)...\n")
cuGetContext(TRUE)
cat("done. Profiling CUDA code...\n")

cat("Loading module...\n")
m = loadModule("truncnorm.ptx")  
cat("done. Extracting kernel...\n")
k = m$rtruncnorm_kernel



###################################################################
"compute_grid" <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
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
###################################################################

K = 7
###
N = 10^K	
###

if(N < 10^8){
	bg = compute_grid(N)
	grid_dims = bg$grid_dims
	block_dims = bg$block_dims
}else{
    threads_per_block <- 512L
	block_dims <- c(threads_per_block, 1L, 1L)
	grid_d1 <- floor(sqrt(N/threads_per_block))
	grid_d2 <- ceiling(N/(grid_d1*threads_per_block))
	grid_dims <- c(grid_d1, grid_d2, 1L)
}

cat("Grid size:\n")
print(grid_dims)
cat("Block size:\n")
print(block_dims)

nthreads <- prod(grid_dims)*prod(block_dims) 
cat("Total number of threads to launch = ",nthreads,"\n")
if (nthreads < N){
    stop("Grid is not large enough...!")
}

# Change values here as needed

lo_in = rep(-15, N)
hi_in = rep(-10, N)
mu_in = rep(0, N)
sigma_in = rep(1, N)
truncnorm_samples = rep(0,N)
maxtries_in = 2000L

cat("Begin CUDA execution...\n")
cuda_copy_to_time = system.time({
	vals_dev = copyToDevice(truncnorm_samples)
	mu_dev = copyToDevice(mu_in)
	sigma_dev = copyToDevice(sigma_in)
	lo_dev = copyToDevice(lo_in)
	hi_dev = copyToDevice(hi_in)
})

cuda_run_time = system.time({
	.cuda(k, vals_dev, N, mu_dev, sigma_dev, lo_dev, hi_dev, N, N, N, N, maxtries_in, gridDim = grid_dims, blockDim = block_dims)
})
cuda_copy_from_time = system.time({
	truncnorm_samples = copyFromDevice(obj=vals_dev,nels=vals_dev@nels,type="float")
})

print(summary(truncnorm_samples))

cat("Cuda time:\n")
print(cuda_run_time)


cpu_time = system.time({
	truncnorm_samples_2 = rtruncnorm(N, 0,1.5, 2,1)
})

print(summary(truncnorm_samples_2))

cat("CPU time:\n")
print(cpu_time)
	
time_table = data.frame(gpu_copy = cuda_copy_to_time[3],
			   			gpu_run = cuda_run_time[3],
			   			gpu_copyback = cuda_copy_from_time[3],
			   			cpu_total = cpu_time[3])
					  
write.csv(time_table, paste0("time_table_K", K, ".csv"))


pdf("tail_sampling.pdf")
plot(density(truncnorm_samples), main = "Density of GPU Truncated Normal Samples")
dev.off()

# Free memory...
rm(list=ls())

q("no")
