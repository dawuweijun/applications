/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include "../../helper.h"

#ifdef USE_DOUBLE
 	double f(double x)
#else
 	float f(float x)
#endif
 	{
 		return x;
 	}


int main(int argc, char **argv){

	if(flagValueInt(argv, argc, "threads") > 0){
		__cilkrts_set_param("nworkers", flagValueText(argv, argc, "threads"));
	}
	else{
		__cilkrts_set_param("nworkers", "2");
	}

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Trapezoidal Rule - Cilk - Threads: %d - Double Precision\n", __cilkrts_get_nworkers());
		#else
			printf("Trapezoidal Rule - Cilk - Threads: %d - Single Precision\n", __cilkrts_get_nworkers());
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}


	#ifdef USE_DOUBLE
		double a, b, I, h;
	#else		 
		float a, b, I, h;
	#endif
	long long int i, n;

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
    
    if(checkFlag(argv, argc, "debug")){
    	printf("Begin main loop\n");
    }

    I = 0;


    #ifdef USE_DOUBLE
    	CILK_C_REDUCER_OPADD(sum, double, 0);
    #else
    	CILK_C_REDUCER_OPADD(sum, float, 0);
    #endif
    CILK_C_REGISTER_REDUCER(sum);

    timer = crono();

    cilk_for (i = 1; i < n; i++) {
    	#ifdef USE_DOUBLE
    		REDUCER_VIEW(sum) += f(a + (double)i * h);
    	#else
    		REDUCER_VIEW(sum) += f(a + (float)i * h);
    	#endif
	}
	timer = crono() - timer;

	I = (h / 2) * (f(a) + 2 * sum.value + f(b));
	CILK_C_UNREGISTER_REDUCER(sum);

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