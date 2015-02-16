/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "../../helper.h"

#ifdef USE_DOUBLE
 	double f(double x)
#else
 	float f(float x)
#endif
 	{
 		return x;
 	}

// Mutex
pthread_mutex_t mutex;

// Threads
long int num_threads;

// Trapezoidal

#ifdef USE_DOUBLE
	double a, b, I, h;
#else		 
	float a, b, I, h;
#endif
	long long int n;


void *trapezoidal(void *param){
	long int id = (long int) param;
	long long int i, chunk = (n - 2) / (long long int) num_threads;
	#ifdef USE_DOUBLE
		double I_local = 0.0;
	#else
		float I_local  = 0.0;
	#endif

	for (i = (long long int) id * chunk; i < (long long int) (id + 1) * chunk; i++) {
	    	#ifdef USE_DOUBLE
	    		I_local += f(a + (double)i * h);
	    	#else
	    		I_local += f(a + (float)i * h);
	    	#endif
	}

	pthread_mutex_lock(&mutex);
	I += I_local;
	pthread_mutex_unlock(&mutex);

	pthread_exit((void *) 0);
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
			printf("Trapezoidal Rule - Pthreads - Threads: %ld - Double Precision\n", num_threads);
		#else
			printf("Trapezoidal Rule - Pthreads - Threads: %ld - Single Precision\n", num_threads);
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

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

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("[%.2lf %.2lf] ---- n = %lld --- h = %.10lf\n", a, b, n, h);
		#else
			printf("[%.2f %.2f] ---- n = %lld --- h = %.10f\n", a, b, n, h);
		#endif
	}

	//Execution time
	double timer = 0.0;

	//Threads
	long int i;
	pthread_t *thread;
	pthread_attr_t attr;
	int error;
	void *status;
    
    if(checkFlag(argv, argc, "debug")){
    	printf("Begin main loop\n");
    }

	/* Posix Threads */
	thread = (pthread_t *) malloc(num_threads * sizeof(pthread_t));

	pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);	

    pthread_mutex_init(&mutex, NULL);

	    I = 0;

	    timer = crono();

		for(i = 0; i < num_threads; ++i){
			if(checkFlag(argv, argc, "debug"))
				printf("\tCreating thread %ld\n", i);
			
			error = pthread_create(&thread[i], &attr, trapezoidal, (void *)i);
			if(error)
				HANDLE_ERROR();
		}

		if(checkFlag(argv, argc, "debug"))
			printf("\n");

		pthread_attr_destroy(&attr);
		for(i = 0; i < num_threads; ++i){
			error = pthread_join(thread[i], &status);
			if(error)
				HANDLE_ERROR();
			if(checkFlag(argv, argc, "debug"))
				printf("\tCompleted join with thread %ld having a status of %ld\n", i, (long int) status);
		}

		timer = crono() - timer;

		I = (h / 2) * (f(a) + 2 * I + f(b));

	/* Posix Threads */
	pthread_mutex_destroy(&mutex);

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

	pthread_exit((void *) 0);
}