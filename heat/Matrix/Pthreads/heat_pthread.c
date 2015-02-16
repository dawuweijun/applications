/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "../../../helper.h"

#define MAX_TEMP 1.0000f
#define MIN_TEMP 0.0001f
#define SPEED    0.2500f

#ifdef USE_DOUBLE
	double **old, **temp;
#else
	float **old, **temp;
#endif

unsigned int size, iterations;

// Barrier
pthread_barrier_t barrier;

// Threads
long int num_threads;

void *heat(void *param){
	long int id = (long int) param;
	unsigned int chunk = size / num_threads;
	unsigned int time;
	unsigned int x, y;

    #ifdef USE_DOUBLE
    	double t, l, c, r, b;	
    #else
   		float  t, l, c, r, b;
   	#endif

	for(time = 1; time <= iterations; time++){

		for(x = id * chunk; x < (id + 1) * chunk; x++){
			for(y = 0; y < size; y++){
				old[x][y] = temp[x][y];
			}
		}
		pthread_barrier_wait(&barrier);

		for(x = id * chunk; x < (id + 1) * chunk; x++){
			for(y = 0; y < size; y++){

 				if(x == 0)
			    	t = old[x][y];
			    else
			    	t = old[x - 1][y];

			    if(y == 0)
			   		l = old[x][y];
			   	else
			   		l = old[x][y - 1];

			   	c = old[x][y];

			   	if(y == size - 1)
			   		r = old[x][y];
			   	else
			   		r = old[x][y + 1];

			   	if(x == size - 1)
			   		b = old[x][y];
			   	else
			   		b = old[x + 1][y];

			    temp[x][y] = c + SPEED * (t + l + r + b - 4 * c);
			}
		}
		pthread_barrier_wait(&barrier);
	}

	pthread_exit((void *) 0);
}

#ifdef USE_DOUBLE
	void print(double **temp, unsigned int size)
#else
	void print(float  **temp, unsigned int size)
#endif
	{
		unsigned int x, y;

		for(x = 0; x < size; x++){
			for(y = 0; y < size; y++){
	    		#ifdef USE_DOUBLE
	    			printf("%.4lf\t", temp[x][y]);
	    		#else
	    			printf("%.4f\t",  temp[x][y]);
	    		#endif
    		}

    		printf("\n");
	    }
	    printf("\n\n");
	}

int main(int argc, char **argv){
	
	if(flagValueInt(argv, argc, "threads") > 0){
		num_threads = flagValueInt(argv, argc, "threads");
	}
	else{
		num_threads = 2;
	}

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Heat Equation - Matrix - Pthreads - Threads: %ld - Double Precision\n", num_threads);
		#else
			printf("Heat Equation - Matrix - Pthreads - Threads: %ld - Single Precision\n", num_threads);
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

	size = 8;
	iterations = 10;

	if(flagValueInt(argv, argc, "size") > 0)
		size = flagValueInt(argv, argc, "size");

	if(flagValueInt(argv, argc, "iterations") > 0)
		iterations = flagValueInt(argv, argc, "iterations");

	if(checkFlag(argv, argc, "debug")){
		printf("Size = %d x %d ---- Iterations = %d\n", size, size, iterations);
	}

	//Threads
	long int i;
	pthread_t *thread;
	pthread_attr_t attr;
	int error;
	void *status;

	unsigned int x, y;

	//Execution time
	double timer = {0.0};

	#ifdef USE_DOUBLE
		old  = (double **) malloc(size * sizeof(double *));
		for(i = 0; i < size; i++)
			old[i]  = (double *) malloc(size * sizeof(double));

		temp = (double **) malloc(size * sizeof(double *));
		for(i = 0; i < size; i++)
			temp[i] = (double *) malloc(size * sizeof(double));
	#else
		old  = (float **) malloc(size * sizeof(float *));
		for(i = 0; i < size; i++)
			old[i]  = (float *) malloc(size * sizeof(float));

		temp = (float **) malloc(size * sizeof(float *));
		for(i = 0; i < size; i++)
			temp = (float *) malloc(size * sizeof(float));
	#endif

    for(x = 0; x < size; x++){
    	for(y = 0; y < size; y++){
    		temp[x][y] = 0;
    	}
    }

    for(x = 0; x < size; x++)
    	for(y = 0; y < size; y++)
    		if(x <= (unsigned int)(0.05 * size) || x >= size - (unsigned int)(0.05 * size) || y <= (unsigned int)(0.05 * size) || y >= size - (unsigned int)(0.05 * size))
    			temp[x][y] = (MAX_TEMP + MIN_TEMP) / 2;
    
    if(checkFlag(argv, argc, "debug")){
    	printf("Begin main loop\n");
    }

/* Posix Threads */
	thread = (pthread_t *) malloc(num_threads * sizeof(pthread_t));

	pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);	

    pthread_barrier_init(&barrier, NULL, num_threads);	

	timer = crono();

	for(i = 0; i < num_threads; ++i){
		if(checkFlag(argv, argc, "debug"))
			printf("Creating thread %ld\n", i);
		
		error = pthread_create(&thread[i], &attr, heat, (void *)i);
		if(error)
			HANDLE_ERROR();
	}

	pthread_attr_destroy(&attr);
	for(i = 0; i < num_threads; ++i){
		error = pthread_join(thread[i], &status);
		if(error)
			HANDLE_ERROR();
		if(checkFlag(argv, argc, "debug"))
			printf("Completed join with thread %ld having a status of %ld\n", i, (long int) status);
	}

	pthread_barrier_destroy(&barrier);

	timer = crono() - timer;

	if(checkFlag(argv, argc, "debug")){
	    printf("Finish main loop\n");
	    printf("Time\n");

	    printf("\texecution time: %.10lf segs.\n", timer);

	    printf("End of the execution\n\n");
	}	

	if(checkFlag(argv, argc, "check") && size == 8 && iterations == 10){
		print(temp, size);

		printf("OK\n\n0.2809	0.2621	0.2322	0.2048	0.1846	0.1744	0.1699	0.1686\n0.2621	0.2407	0.2108	0.1791	0.1590	0.1467	0.1422	0.1406\n0.2322	0.2108	0.1743	0.1426	0.1174	0.1052	0.0993	0.0977\n0.2048	0.1791	0.1426	0.1043	0.0791	0.0641	0.0582	0.0561\n0.1846	0.1590	0.1174	0.0791	0.0504	0.0354	0.0286	0.0266\n0.1744	0.1467	0.1052	0.0641	0.0354	0.0193	0.0126	0.0104\n0.1699	0.1422	0.0993	0.0582	0.0286	0.0126	0.0057	0.0036\n0.1686	0.1406	0.0977	0.0561	0.0266	0.0104	0.0036	0.0015\n");
	}
	
	for(i = 0; i < size; i++)
		free(temp[i]);
    free(temp);

    for(i = 0; i < size; i++)
		free(old[i]);
    free(old);

    fprintf(stderr, "%.10lf\n", timer);

	return 0;
}