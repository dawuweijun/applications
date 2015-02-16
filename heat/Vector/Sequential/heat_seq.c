/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include "../../../helper.h"

#define MAX_TEMP 1.0000f
#define MIN_TEMP 0.0001f
#define SPEED    0.2500f

enum Heat{COPY = 0, BLEND, TEMP, TOTAL};

#ifdef USE_DOUBLE
	void copy(double *temp, double *old, unsigned int size)
#else
	void copy(float  *temp, float  *old, unsigned int size)
#endif
	{
		unsigned int offset;

		for(offset = 0; offset < size * size; offset++){
			old[offset] = temp[offset];
		}

	}

#ifdef USE_DOUBLE
	void blend(double *temp, double *old, unsigned int size)
#else
	void blend(float  *temp, float  *old, unsigned int size)
#endif
	{
	    unsigned int x, y, offset;
	    unsigned int left, right, top, bottom;

	    #ifdef USE_DOUBLE
	    	double t, l, c, r, b;	
	    #else
	   		float  t, l, c, r, b;
	   	#endif
	    
		for(offset = 0; offset < size * size; offset++){
			x = offset / size;
			y = offset % size;

			left  = offset - 1;
		    right = offset + 1;
		    
		    if (y == 0)
		    	left++;
		    if (y == size - 1)
		    	right--; 

		    top    = offset - size;
		    bottom = offset + size;
		    
		    if (x == 0)
		    	top    += size;
		    if (x == size - 1)
		    	bottom -= size;

		   	t = old[top];
		   	l = old[left];
		   	c = old[offset];
		   	r = old[right];
		   	b = old[bottom];	

		    temp[offset] = c + SPEED * (t + l + r + b - 4 * c);
	    }
	}

#ifdef USE_DOUBLE
	void print(double *temp, unsigned int size)
#else
	void print(float  *temp, unsigned int size)
#endif
	{
		unsigned int offset;

		for(offset = 0; offset < size * size; offset++){
    		#ifdef USE_DOUBLE
    			printf("%.4lf\t", temp[offset]);
    		#else
    			printf("%.4f\t",  temp[offset]);
    		#endif

    			if((offset + 1) % size == 0 && offset > 0)
    				printf("\n");
	    }
	    printf("\n\n");
	}

int main(int argc, char **argv){

	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Heat Equation - Vector - Sequential - Double Precision\n");
		#else
			printf("Heat Equation - Vector - Sequential - Single Precision\n");
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

	unsigned int size = 8;
	unsigned int iterations = 10;

	if(flagValueInt(argv, argc, "size") > 0)
		size = flagValueInt(argv, argc, "size");

	if(flagValueInt(argv, argc, "iterations") > 0)
		iterations = flagValueInt(argv, argc, "iterations");

	if(checkFlag(argv, argc, "debug")){
		printf("Size = %d x %d ---- Iterations = %d\n", size, size, iterations);
	}

	unsigned int offset, i, x, y;

	//Execution time
	double timer[4] = {0.0};

	#ifdef USE_DOUBLE
		double *old  = (double *) malloc(size * size * sizeof(double));
		double *temp = (double *) malloc(size * size * sizeof(double));
	#else
		float *old  = (float *)   malloc(size * size * sizeof(float));
		float *temp = (float *)   malloc(size * size * sizeof(float));
	#endif

    for (offset = 0; offset < size * size; offset++){
    	temp[offset] = 0;
    }

    for(x = 0; x < size; x++)
    	for(y = 0; y < size; y++)
    		if(x <= (unsigned int)(0.05 * size) || x >= size - (unsigned int)(0.05 * size) || y <= (unsigned int)(0.05 * size) || y >= size - (unsigned int)(0.05 * size))
    			temp[x * size + y] = (MAX_TEMP + MIN_TEMP) / 2;
    
    if(checkFlag(argv, argc, "debug")){
    	printf("Begin main loop\n");
    }

    timer[TOTAL] = crono();
    for (i = 1; i <= iterations; i++) {
    	timer[TEMP] = crono();
		copy(temp, old, size);
		timer[TEMP] = crono() - timer[TEMP];
		timer[COPY] += timer[TEMP];

		timer[TEMP] = crono();
		blend(temp, old, size);
		timer[TEMP] = crono() - timer[TEMP];
		timer[BLEND] += timer[TEMP ];
	}
	timer[TOTAL] = crono() - timer[TOTAL];

	if(checkFlag(argv, argc, "debug")){
	    printf("Finish main loop\n");
	    printf("Time\n");

	    printf("\tcopy time: %.10lf segs.\n", timer[COPY]);

	    printf("\tblend time: %.10lf segs.\n", timer[BLEND]);

	    printf("\texecution time: %.10lf segs.\n", timer[TOTAL]);

	    printf("End of the execution\n\n");
	}	

	if(checkFlag(argv, argc, "check") && size == 8 && iterations == 10){
		print(temp, size);

		printf("OK\n\n0.2809	0.2621	0.2322	0.2048	0.1846	0.1744	0.1699	0.1686\n0.2621	0.2407	0.2108	0.1791	0.1590	0.1467	0.1422	0.1406\n0.2322	0.2108	0.1743	0.1426	0.1174	0.1052	0.0993	0.0977\n0.2048	0.1791	0.1426	0.1043	0.0791	0.0641	0.0582	0.0561\n0.1846	0.1590	0.1174	0.0791	0.0504	0.0354	0.0286	0.0266\n0.1744	0.1467	0.1052	0.0641	0.0354	0.0193	0.0126	0.0104\n0.1699	0.1422	0.0993	0.0582	0.0286	0.0126	0.0057	0.0036\n0.1686	0.1406	0.0977	0.0561	0.0266	0.0104	0.0036	0.0015\n");
	}
	
    free(temp);
    free(old);

    fprintf(stderr, "%.10lf %.10lf %.10lf\n", timer[COPY], timer[BLEND], timer[TOTAL]);

	return 0;
}