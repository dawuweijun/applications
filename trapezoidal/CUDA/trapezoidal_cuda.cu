/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include "../../boltzmann/CUDA/helper.cuh"
#include <cuda_runtime.h>

static void GPUHandleError( cudaError_t err, const char *file, const int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

#define GPU_HANDLE_ERROR( err ) (GPUHandleError( err, __FILE__, __LINE__ ))

	long long int n;
	int threadsPerBlock, blocksPerGrid;

#ifdef USE_DOUBLE
 	__host__ __device__ double f(double x)
#else
 	__host__ __device__ float f(float x)
#endif
 	{
 		return x;
 	}


#ifdef USE_DOUBLE
	__global__ void sum_reduct(double a, double h, double *sum, long long int n)
#else
	__global__ void sum_reduct(float a, float h, float *sum, long long int n)
#endif
	{
		#ifdef USE_DOUBLE
			__shared__ double cache[256];
			double temp = 0.0;
		#else
			__shared__ float cache[256];
			float temp = 0.0;
		#endif

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int cacheIdx = threadIdx.x;

		while(tid < n){
			temp += f(a + tid * h);
			tid += blockDim.x * gridDim.x;
		}

		// set the cache values
		cache[cacheIdx] = temp;

		// synchronize threads in this block
		__syncthreads();

		//for reductions, threadsPerBlock must be a power of 2
		int i = blockDim.x / 2;
		while(i != 0 ){
			if(cacheIdx < i)
				cache[cacheIdx] += cache[cacheIdx + i];
			__syncthreads();
			i /= 2;
		}

		if(cacheIdx == 0)
			sum[blockIdx.x] = cache[0];

	}

int main(int argc, char **argv){

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Trapezoidal Rule - CUDA - Double Precision\n");
		#else
			printf("Trapezoidal Rule - CUDA - Single Precision\n");
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}


	#ifdef USE_DOUBLE
		double a, b, I, h;
	#else		 
		float a, b, I, h;
	#endif

	if(flagValueLong(argv, argc, "n") > 0)
		n = flagValueLong(argv, argc, "n");
	else
		n = 10;

	if(flagValueReal(argv, argc, "a") > 0){
		a = flagValueReal(argv, argc, "a");
	}
	else
		a = 1;

	if(flagValueReal(argv, argc, "b") > 0)
		b = flagValueReal(argv, argc, "b");
	else
		b = 2;

	h = (b - a) / n;

	//Execution time
	float timer = 0.0;
    
    if(checkFlag(argv, argc, "debug")){
    	printf("Begin main loop\n");
    }

   	// Workstation UNIPAMPA: 0 to Tesla C2075, 1 to Quadro 5000
	GPU_HANDLE_ERROR(cudaSetDevice(0));

	// Device memory and Host memory
	#ifdef USE_DOUBLE
		double *h_sum, *d_sum;
	#else
		float  *h_sum, *d_sum;
	#endif

	// Blocks and Grids
	threadsPerBlock = 256; // remember to chance in cuda kernel
	blocksPerGrid = 32;

	if(((n - 2 + threadsPerBlock - 1) / threadsPerBlock) < 32)
		blocksPerGrid = (n - 2 + threadsPerBlock - 1) / threadsPerBlock;

	if(checkFlag(argv, argc, "debug")){
		printf("\n");
		#ifdef USE_DOUBLE
			printf("[%.2lf %.2lf] ---- n = %lld --- h = %.10lf\n", a, b, n, h);
		#else
			printf("[%.2f %.2f] ---- n = %lld --- h = %.10f\n", a, b, n, h);
		#endif
		printf("threadsPerBlock = %d --- blocksPerGrid = %d\n\n", threadsPerBlock, blocksPerGrid);
	}



	// Timer
	cudaEvent_t start, stop;

	GPU_HANDLE_ERROR(cudaEventCreate(&start));
	GPU_HANDLE_ERROR(cudaEventCreate(&stop));

	// Memory alloc for host and device
	#ifdef USE_DOUBLE	
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_sum,  blocksPerGrid * sizeof(double)));
		h_sum = (double *) malloc(blocksPerGrid * sizeof(double));
	#else
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_sum,  blocksPerGrid * sizeof(float)));
		h_sum = (double *) malloc(blocksPerGrid * sizeof(float));
	#endif

    I = 0;

    // CUDA Code
	GPU_HANDLE_ERROR(cudaDeviceSynchronize());

	GPU_HANDLE_ERROR(cudaEventRecord(start, 0));
		sum_reduct<<<blocksPerGrid, threadsPerBlock>>>(a, h, d_sum, n - 2);
	GPU_HANDLE_ERROR(cudaEventRecord(stop, 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop));

  	GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());

    GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer, start, stop));
    
    // Converting miliseconds to seconds
	timer /= 1000;

	// Stop events
  	GPU_HANDLE_ERROR(cudaEventDestroy(start));
	GPU_HANDLE_ERROR(cudaEventDestroy(stop));

	// Copy device to host
	#ifdef USE_DOUBLE
		GPU_HANDLE_ERROR(cudaMemcpy(h_sum, d_sum, blocksPerGrid * sizeof(double), cudaMemcpyDeviceToHost));
	#else
		GPU_HANDLE_ERROR(cudaMemcpy(h_sum, d_sum, blocksPerGrid * sizeof(float),  cudaMemcpyDeviceToHost));
	#endif

	// Synchronize
	GPU_HANDLE_ERROR(cudaDeviceSynchronize());

	// finish up on the CPU side
	int i;
	I = 0;
	for(i = 0; i < blocksPerGrid; i++)
		I += h_sum[i];

	// CUDA Code again
	// Free device memory
	GPU_HANDLE_ERROR(cudaFree(d_sum));

	// Free host memory
	free(h_sum);

	// Device Reset
	GPU_HANDLE_ERROR(cudaDeviceReset());

	if(checkFlag(argv, argc, "debug"))
		#ifdef USE_DOUBLE
			printf("GPU_sum = %.10lf\n", I);
		#else
			printf("GPU_sum = %.10f\n", I);
		#endif

	I = (h / 2) * (f(a) + 2 * I + f(b));

	if(checkFlag(argv, argc, "debug")){
	    printf("Finish main loop\n");

	    #ifdef USE_DOUBLE
	    	printf("\nTrapezoidal = %.10lf\n", I);
	    #else
	    	printf("\nTrapezoidal = %.10f\n", I);
	    #endif

	    printf("Time\n");

	    printf("\texecution time: %.10lf segs.\n", timer);

	    printf("End of the execution\n\n");
	}	

	if(checkFlag(argv, argc, "check") && a == 1 && b == 2){
		#ifdef USE_DOUBLE
	    	printf("%.10lf\n", I);
	    #else
	    	printf("%.10f\n", I);
	    #endif
		printf("OK\n1.5\n\n");
	}

    fprintf(stderr, "%.10lf\n", timer);

	return 0;
}