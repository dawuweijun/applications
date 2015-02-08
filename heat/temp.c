/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>

#define MAX_TEMP 1.0f
#define MIN_TEMP 0.0001f
#define SPEED   0.25f
#define DIM 10000
#define ITERATIONS 10

void copy(double *temp, double *old){
	int i;
	#pragma omp parallel for default(shared) private(i)
	for(i = 0; i < DIM * DIM; i++)
		old[i] = temp[i];
}

void blend(double *temp, double *old){
    int x, y;
    #pragma omp parallel for default(shared) private(x, y)
    for(x = 0; x < DIM; x++){
    	for(y = 0; y < DIM; y++){
    		int offset = x * DIM + y;
    		int left = offset - 1;
		    int right = offset + 1;
		    if (y == 0)   left++;
		    if (y == DIM-1) right--; 

		    int top = offset - DIM;
		    int bottom = offset + DIM;
		    if (x == 0)   top += DIM;
		    if (x == DIM-1) bottom -= DIM;

		    float t = old[top];
		    float l = old[left];
		    float c = old[offset];
		    float r = old[right];
		    float b = old[bottom];

		    temp[offset] = c + SPEED * (t + b + r + l - 4 * c);
    	}
    }
}

void print(double *temp){
	int x, y;
	for(x = 0; x < DIM; x++){
    	for(y = 0; y < DIM; y++){
    		printf("%.4f\t", temp[x*DIM + y]);
    	}
    	printf("\n");
    }	
    printf("\n\n");
}

int main(){
	int i, x, y;

	double *old = (double*)malloc(DIM * DIM * sizeof(double));
	double *temp = (double*)malloc(DIM * DIM * sizeof(double));
    


    for (i=0; i<DIM*DIM; i++)
    	temp[i] = 0;

    for(x = 0; x < DIM; x++)
    	for(y = 0; y < DIM; y++)
    		if(x <= (int)(0.05 * DIM) || x >= DIM - (int)(0.05 * DIM) || y <= (int)(0.05 * DIM) || y >= DIM - (int)(0.05 * DIM))
    			temp[x * DIM + y] = (MAX_TEMP + MIN_TEMP) / 2;
    
    //print(temp);
    
    for (i = 1; i <= ITERATIONS; i++) {
		copy(temp, old);
		blend(temp, old);
	}
	//print(temp);

	//printf("OK\n\n");

	//printf("0.2809	0.2621	0.2322	0.2048	0.1846	0.1744	0.1699	0.1686\n0.2621	0.2407	0.2108	0.1791	0.1590	0.1467	0.1422	0.1406\n0.2322	0.2108	0.1743	0.1426	0.1174	0.1052	0.0993	0.0977\n0.2048	0.1791	0.1426	0.1043	0.0791	0.0641	0.0582	0.0561\n0.1846	0.1590	0.1174	0.0791	0.0504	0.0354	0.0286	0.0266\n0.1744	0.1467	0.1052	0.0641	0.0354	0.0193	0.0126	0.0104\n0.1699	0.1422	0.0993	0.0582	0.0286	0.0126	0.0057	0.0036\n0.1686	0.1406	0.0977	0.0561	0.0266	0.0104	0.0036	0.0015\n");

    free(temp);
    free(old);

	return 0;
}