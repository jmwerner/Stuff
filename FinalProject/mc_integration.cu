#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{
	
	__device__ float normal_pdf(float x){
		float returner = 0;
		returner = expf(-1 * (pow(x,2)) / 2) / sqrt(2 * (3.141592653));
		return returner;
	}// normal_pdf
	
	__device__ float exp_pdf(float x, float lambda){
		float returner = 0;
		returner = lambda * expf(-1 * lambda * x);
		return returner;
	}// exp_pdf
	

	__constant__ int seed_a = 1234, seed_b = 1423, seed_c = 1842;
	
	__global__ void 
	mc_integration_normal_kernel(int *vals, int vals_length, int n, 
								 float quantile, float sigma_out)
	{
	
		// Usual block/thread indexing...
		int myblock = blockIdx.x + blockIdx.y * gridDim.x;
		int blocksize = blockDim.x * blockDim.y * blockDim.z;
		int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
		int idx = myblock * blocksize + subthread;
	
		// Declaration of other needed variables
		float x_uni = 0, y_uni = 0; 
	
		// Only evaluate if thread is within the range of values we need
		if(idx < n){
			// Setup the RNG:
			curandState rng;
			curand_init(seed_a + idx * seed_b, seed_c, 0, &rng);

			// Provides x uniform on quantile to sigma_out
		   	 x_uni = (sigma_out - quantile) * curand_uniform(&rng) + quantile; 
		    // Provides y uniform on 0 to 1
			y_uni = curand_uniform(&rng);
			
			if(y_uni < normal_pdf(x_uni)){
				vals[idx] = 1;
			}else{
				vals[idx] = 0;
			}//ifelse
		}// if(idx < n)
		return;
	} // mc_integration_normal_kernel
	
	
	__global__ void 
	mc_integration_normal_vegas_kernel(float *vals, int vals_length, int n, 
								 float quantile, float lambda)
	{
		// Usual block/thread indexing...
		int myblock = blockIdx.x + blockIdx.y * gridDim.x;
		int blocksize = blockDim.x * blockDim.y * blockDim.z;
		int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
		int idx = myblock * blocksize + subthread;
	
		// Declaration of other needed variables
		float x_exp = 0, testval = 0;
	
		// Only evaluate if thread is within the range of values we need
		if(idx < n){
			// Setup the RNG:
			curandState rng;
			curand_init(seed_a + idx * seed_b, seed_c, 0, &rng);

			// Provides random value from exponential(lambda)
		   	 x_exp = -log(curand_uniform(&rng))/lambda;

			
			testval = normal_pdf(x_exp) / exp_pdf(x_exp, lambda);
			
			if(x_exp > quantile){
				vals[idx] = testval;
			}else{
				vals[idx] = 0;
			}//ifelse
		}// if(idx < n)
		return;
	} // mc_integration_normal_vegas_kernel
	
	__global__ void 
	mc_integration_normal_kernel_2(float *vals, int vals_length, int N, 
								 float *quantiles, int quantile_num, float sigma_out)
	{
	
		// Usual block/thread indexing...
		int myblock = blockIdx.x + blockIdx.y * gridDim.x;
		int blocksize = blockDim.x * blockDim.y * blockDim.z;
		int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
		int idx = myblock * blocksize + subthread;
	
		// Declaration of other needed variables
		float x_uni = 0, y_uni = 0, quantile; 
		int counter = 0, i = 0;
	
		// Only evaluate if thread is within the range of values we need
		if(idx < quantile_num){
			// Setup the RNG:
			curandState rng;
			curand_init(seed_a + idx * seed_b, seed_c, 0, &rng);
		
			quantile = quantiles[idx];
			
			for(i = 0; i < N; i++){	

				// Provides x uniform on quantile to sigma_out
		   		x_uni = (sigma_out - quantile) * curand_uniform(&rng) + quantile; 
			    // Provides y uniform on 0 to 1
				y_uni = curand_uniform(&rng);
			
				if(y_uni < normal_pdf(x_uni)){
					counter++;
				}// if
			}// for
			
			vals[idx] = counter * (sigma_out - quantile) / N;
			
		}// if(idx < n)
		return;
	} // mc_integration_normal_kernel
	
	
} // END extern "C"


