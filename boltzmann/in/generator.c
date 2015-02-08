/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
	int i, j, k, max_x, max_y, max_z;

	if(argc == 2){
		max_x = atoi(argv[1]);
		max_y = max_x;
		max_z = max_x;
	}
	else if(argc == 4){
		max_x = atoi(argv[1]);
		max_y = atoi(argv[2]);
		max_z = atoi(argv[3]);
	}else{
		fprintf(stderr, "Usage either:\n\t%s [size_x] [size_y] [size_z] > output.obs\n", argv[0]);
		fprintf(stderr, "or:\n\t%s [size] > output.obs\n\n", argv[0]);
		return 0;
	}

	printf("%d %d %d\n", max_x, max_y, max_z);
	
	for(i = 1; i < max_x; i++)
		for(j = 1; j < max_y; j++)
			printf("%d %d %d\n", i, j, 1);

	for(i = 1; i < max_x; i++)
		for(j = 1; j < max_y; j++)
			printf("%d %d %d\n", i, j, max_z - 1);

	for(i = 1; i < max_x; i++)
		for(k = 1; k < max_z; k++)
			printf("%d %d %d\n", i, 1, k);

	for(i = 1; i < max_x; i++)
		for(k = 1; k < max_z; k++)
			printf("%d %d %d\n", i, max_y - 1, k);

	for(j = 1; j < max_y; j++)
		for(k = max_z/3; k < 2*max_z/3; k++)
			printf("%d %d %d\n", max_x/3, j, k);

	return 0;
}