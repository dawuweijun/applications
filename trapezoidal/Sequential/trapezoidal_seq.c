/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
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

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Trapezoidal Rule - Sequential - Double Precision\n");
		#else
			printf("Trapezoidal Rule - Sequential - Single Precision\n");
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

    timer = crono();
    for (i = 1; i < n; i++) {
    	#ifdef USE_DOUBLE
    		I += f(a + (double)i * h);
    	#else
    		I += f(a + (float)i * h);
    	#endif
	}
	timer = crono() - timer;

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