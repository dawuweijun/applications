/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "lb_3D.h"
#include "helper.h"

/* Alloc memory space to the grid */
void alloc_lattice(s_lattice *l, unsigned short int debug){
	
	l->obst = (unsigned short int *) calloc(l->lx * l->ly * l->lz, sizeof(unsigned short int **));
	#ifdef USE_DOUBLE
		l->node = (double *) calloc(l->lx * l->ly * l->lz * NDIM , sizeof(double***));
		l->temp = (double *) calloc(l->lx * l->ly * l->lz * NDIM , sizeof(double***));
	#else
		l->node = (float *) calloc(l->lx * l->ly * l->lz * NDIM , sizeof(float***));
		l->temp = (float *) calloc(l->lx * l->ly * l->lz * NDIM , sizeof(float***));		
	#endif

	if(debug)
		printf("Alloc lattices sucessful\n");
}

/* Returns velocity of middle element in axis x*/
#ifdef USE_DOUBLE
	double calc_velocity(s_lattice *l, unsigned short int debug)
#else
	float calc_velocity(s_lattice *l, unsigned short int debug)
#endif
	{
		unsigned short int x, y, z, i, n_free;
		#ifdef USE_DOUBLE
			double u_x, d_loc, velocity;
		#else
			float u_x, d_loc, velocity;
		#endif

		x = l->lx / 2 + l->lx % 2; /* Middle element */
		n_free = 0;
		u_x = 0;

		for(y = 0; y < l->ly; y++)
			for(z = 0; z < l->lz; z++)
				if (l->obst[x * l->ly * l->lz + y * l->lz + z] == false){
					d_loc = 0;
					for (i = 0; i < NDIM; i++)
						d_loc = d_loc + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + i];
					u_x += (l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14]  - (l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13])) / d_loc;
					n_free++;
				}

		velocity = u_x / n_free;
		if(debug)
			#ifdef USE_DOUBLE
				printf("\tvelocity: %.10lf -0.1528171488 \n", velocity);
			#else
				printf("\tvelocity: %.10f -0.1528169662 \n", velocity);
			#endif

		if(velocity + 0.136921 > 0.1){
			printf("MOTHERFUCK ERROR!\n");
			HANDLE_ERROR();
		}
		
		return velocity;
	}

/* Returns sum of densities in all nodes */
double check_density(s_lattice *l, unsigned short int debug){
	unsigned short int x, y, z, n;
	#ifdef USE_DOUBLE
		double sum = 0.0;
	#else
		float sum = 0.0;
	#endif

	for (x = 0; x < l->lx; x++)
		for (y = 0; y < l->ly; y++)
			for (z = 0; z < l->lz; z++)
				for (n = 0; n < NDIM; n++)
					sum += l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + n];

	if(debug)
		#ifdef USE_DOUBLE
			printf("\tIntegral density = %.10lf\n", sum);
		#else
			printf("\tIntegral density = %.10f\n", sum);
		#endif

	return sum;
}

/* Returns reynolds number of simulation in specific iteration(time) */
#ifdef USE_DOUBLE
	double comp_rey(s_lattice *l, int time, unsigned short int debug)
#else
	float comp_rey(s_lattice *l, int time, unsigned short int debug)
#endif
	{
		#ifdef USE_DOUBLE
			double vel, visc, rey;
		#else
			float vel, visc, rey;
		#endif
		
		vel = calc_velocity(l, debug);
		visc = 1.0 / 6.0 * (2.0 / OMEGA - 1.0);
		rey = vel * REYNOLDS / visc;

		if(debug){
			printf("Compute Reynolds number sucessful\n");
			printf("\titerations: %d\n", time);
			#ifdef USE_DOUBLE
				printf("\tviscosity: %.10lf\n", visc);
				printf("\taverage velocity: %.10lf\n", vel);
				printf("\treynolds number: %.10lf\n", rey);
			#else
				printf("\tviscosity: %.10f\n", visc);
				printf("\taverage velocity: %.10f\n", vel);
				printf("\treynolds number: %.10f\n", rey);				
			#endif
		}

		return rey;
	}

/* Free memory space of the grid */
void dealloc_lattice(s_lattice *l, unsigned short int debug){
	
	free(l->obst);
	free(l->node);
	free(l->temp);
	
	if(debug)
		printf("Dealloc lattices sucessful\n");
}

/* Initializes lattice node with three density levels */
void init_density(s_lattice * l, unsigned short int debug){
	//local variables
	unsigned short int x, y, z, n;
	#ifdef USE_DOUBLE
		double C[NDIM];
	#else
		float C[NDIM];
	#endif
	
	//Coefficients for the particles in XY-plane
	C[0] = DENSITY * T0;
	C[1] = DENSITY * T1;
	C[2] = DENSITY * T2;
	C[3] = DENSITY * T1;
	C[4] = DENSITY * T2;
	C[5] = DENSITY * T1;
	C[6] = DENSITY * T2;
	C[7] = DENSITY * T1;
	C[8] = DENSITY * T2;

	//Coefficients for the particles in XZ-plane
	C[9]  = DENSITY * T2;
	C[10] = DENSITY * T1;
	C[11] = DENSITY * T2;
	C[12] = DENSITY * T2;
	C[13] = DENSITY * T1;
	C[14] = DENSITY * T2;

	//Coefficients for the particles in YZ-plane
	C[15] = DENSITY * T2;
	C[16] = DENSITY * T2;
	C[17] = DENSITY * T2;
	C[18] = DENSITY * T2;

	for (x = 0; x < l->lx; x++)
		for (y = 0; y < l->ly; y++)
			for (z = 0; z < l->lz; z++)
				for (n = 0; n < NDIM; n++)
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + n] = C[n];

	if(debug)
		printf("Init density sucessful\n");

}

/* Print parameters from properties structure */
void print_parameters(char **argv, int argc){	
	if(checkFlag(argv, argc, "debug")){
		printf("Parameters read sucessful\n");
		printf("\tnumber of iterations: %d\n", ITERATIONS);
		#ifdef USE_DOUBLE
			printf("\tfluid density per link: %.10lf\n", DENSITY);
			printf("\tdensity redistribution: %.10lf\n", ACCELERATION);
			printf("\trelaxation parameter: %.10lf\n", OMEGA);
		#else
			printf("\tfluid density per link: %.10f\n", DENSITY);
			printf("\tdensity redistribution: %.10f\n", ACCELERATION);
			printf("\trelaxation parameter: %.10f\n", OMEGA);
		#endif
		printf("\tlinear dimension: %d.0\n", REYNOLDS);
	}
}

/* Read obstacles from file and save in lattice structure */
s_lattice *read_obstacles(char **argv, int argc){ 
	unsigned short int i, j, k, max = 0;
	s_lattice *l = (s_lattice *) malloc(sizeof(s_lattice));
	
	if(!l){
		fprintf(stderr, "Could not alloc s_lattice\n");
		HANDLE_ERROR();
	}

	FILE *archive = fopen(flagValueText(argv, argc, "obstacle"), "r");
	if(!archive){
		fprintf(stderr, "Could not read colision input file\nFile: %s\n\n", flagValueText(argv, argc, "obstacle"));
		HANDLE_ERROR();
	}

	if(!fscanf(archive, "%hu", &l->lx)){
		fprintf(stderr, "Could not read nodes number in axis x\n");
		HANDLE_ERROR();
	}

	if(!fscanf(archive, "%hu", &l->ly)){
		fprintf(stderr, "Could not read nodes number in axis y\n");
		HANDLE_ERROR();	
	}

	if(!fscanf(archive, "%hu", &l->lz)){
		fprintf(stderr, "Could not read nodes number in axis z\n");
		HANDLE_ERROR();	
	}
	
	alloc_lattice(l, checkFlag(argv, argc, "debug"));

	/* Reading obstacle points */
	while(fscanf(archive, "%hu %hu %hu", &i, &j, &k) != EOF){
		//Check if i, j and k are less then x_max, y max and z_max
		if(i < l->lx && j < l->ly && j < l->lz)
			l->obst[i * l->ly * l->lz + j * l->lz + k] = true;
		else{
			fprintf(stderr, "Obstacle point[%hu,%hu,%hu] invalid!\n\n", i, j, k);
		  	HANDLE_ERROR();	
		}
		++max;
	}
	l->n_obst = max;
	fclose(archive);

	if(checkFlag(argv, argc, "debug")){
		printf("Obstacles read sucessful\n");
		printf("\tnodes number in axis x: %u\n", l->lx);
		printf("\tnodes number in axis y: %u\n", l->ly);
		printf("\tnodes number in axis z: %u\n", l->lz);
		printf("\tlattice dimension elements: %u\n", NDIM);
		printf("\tobstacles amount: %u\n", max);
	}
	
	return l;
}

/* Write velocities in x and y, and press in a file. */
/* In the future, if matlab is installed, save eps figures of simulation. */
void write_results(char *file, s_lattice *l){
		unsigned short int x, y, z, i;
		unsigned short int obsval;
		#ifdef USE_DOUBLE
			double u_x, u_y, u_z, d_loc, press;	
			double c_squ = 1.0 / 3.0; /* Square speed of sound */
		#else
			float u_x, u_y, u_z, d_loc, press;	
			float c_squ = 1.0 / 3.0; /* Square speed of sound */
		#endif

		FILE *archive = fopen(file, "w");
		if(!archive){
			fprintf(stderr, "Could not write output file.\nFile: %s\n\n", file);
			HANDLE_ERROR();
		}

		fprintf(archive,"VARIABLES = X, Y, Z, VX, VY, VZ, PRESS, OBST\n");
		fprintf(archive,"ZONE I= %d, J= %d, K= %d F=POINT\n", l->lx, l->ly, l->lz);

		for(z = 0; z < l->lz; z++){
			for(y = 0; y < l->ly; y++){
				for(x = 0; x < l->lx; x++){
					if (l->obst[x * l->ly * l->lz + y * l->lz + z] == true){
						u_x = 0.0;
						u_y = 0.0;
						u_z = 0.0;

						obsval = true;
						press = DENSITY * c_squ; /* average pressure */
					} else {
						d_loc = 0.0;
						for(i = 0; i < NDIM; i++)
							d_loc += l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + i]; /* integral local density */
						
						//x-, y- and z- velocity components
						u_x = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14]  - (l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13]);
						u_x /= d_loc;

						u_y = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18]  - (l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17]);
						u_y /= d_loc;
						
						u_z = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16]  - (l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] + l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18]);
						u_z /= d_loc;

						obsval = false;
						press = d_loc * c_squ;
					}
					#ifdef USE_DOUBLE
						fprintf(archive,"%d %d %d %.10lf %.10lf %.10lf %.10lf %d\n", x, y, z, u_x, u_y, u_z, press, obsval);
					#else
						fprintf(archive,"%d %d %d %.10lf %.10f %.10f %.10f %d\n", x, y, z, u_x, u_y, u_z, press, obsval);
					#endif
				}
			}
		}
		fclose(archive);
	}