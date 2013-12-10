//Kernel for computing truncated normal samples

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{
	
	// Based on example code for random exponential
	__device__ float rexpo(curandState *state, float lambda){
		float value;
		value = -log(curand_uniform(state))/lambda;
		return value;
	}// rexpo
	
	__device__ float psi_function(float lo, float alph, float z){
		float returner = 0;
		if(lo < alph){
		returner = expf(-1*(pow((alph - z),2)/2));
		}else{
		returner = expf(-1*pow((lo-alph),2)/2)*expf(-1*(pow((alph - z),2)/2));
		}
		return returner;
	}// psi_function

	__constant__ int seed_a = 1234, seed_b = 1423, seed_c = 1842;
	
	__global__ void 
	rtruncnorm_kernel(float *vals, int n, 
			  float *mu, float *sigma, 
			  float *lo, float *hi,
        		  int mu_len, int sigma_len,
      			  int lo_len, int hi_len,
	 		  int maxtries)
	{
	
		// Usual block/thread indexing...
		int myblock = blockIdx.x + blockIdx.y * gridDim.x;
		int blocksize = blockDim.x * blockDim.y * blockDim.z;
		int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
		int idx = myblock * blocksize + subthread;
	
		// Declaration of other needed variables
		int found = 0, right_truncation = 0, while_iter = 0;
		float testval = 0, lowbound = 0, zval = 0, alpha = 0, psi = 0, u = 0;
	
		// Only evaluate if thread is within the range of values we need
		if(idx < n){
			// Setup the RNG:
			curandState rng;
			curand_init(seed_a + idx * seed_b,seed_c,0,&rng);
			
	    	// Sample:
		while((found == 0) && while_iter < maxtries){
			testval = mu[idx] + sigma[idx] * curand_normal(&rng);
			if(testval < hi[idx]){
				if(testval > lo[idx]){
				found = 1;
			    }
			}
			while_iter++;
		}//while
		
		if(found == 1){
			vals[idx] = testval;
		}else{
			// Robert's approximation method - coded for right truncation only
			if(lo[idx] > mu[idx]){
			right_truncation = 1;
			lowbound = lo[idx];
			}else{
			lowbound = -1 * hi[idx] + 2 * mu[idx];
			}
	
			// Compute optimum alpha (does not change and doesn't need to be computed in loop below)
	
			alpha = (lowbound + sqrtf(pow(lowbound,2)+4))/2;
			while(!found){
			    zval = lowbound + rexpo(&rng, alpha);
			    psi = psi_function(lowbound, alpha, zval);
			    u = curand_uniform(&rng);
			    if(u < psi){
				found = 1;
				if(right_truncation == 1){
				    vals[idx] = mu[idx] + sigma[idx] * zval;
				}else{
			    	    vals[idx] = mu[idx] - sigma[idx] * zval;
				}
			    }		
			}// while
		}//ifelse
	
			
		}// if(idx < n)
		return;
	} // rtruncnorm_kernel

} // END extern "C"


