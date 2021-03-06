/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "../../../helper.h"
#include "../D3Q19/lb_3D.h"

// Lattice structure
s_lattice *lattice;

// Barrier
pthread_barrier_t barrier;

// Threads
long int num_threads;
unsigned int *chunks;

void *boltzmann(void *param){
	s_lattice *l = lattice;
	int time;
	unsigned short int x, y, z;
	unsigned short int x_e, x_w, y_u, y_d, z_n, z_s;
	#ifdef USE_DOUBLE
		double u_x, u_y, u_z;
		double u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
	#else
		float u_x, u_y, u_z;
		float u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
	#endif
	long int id = (long int) param;

	for(time = 1; time <= ITERATIONS; time++){

		for(x = chunks[id]; x < chunks[id + 1]; ++x)
			for(y = 0; y < l->ly; y++)
				for(z = 0; z < l->lz; z++)
					if(l->obst[x][y][z] == false){
						l->node[x][y][z][1]  += T_1;
						l->node[x][y][z][2]  += T_2;
						l->node[x][y][z][4]  -= T_2;
						l->node[x][y][z][5]  -= T_1;
						l->node[x][y][z][6]  -= T_2;
						l->node[x][y][z][8]  += T_2;
						l->node[x][y][z][9]  += T_2;
						l->node[x][y][z][11] -= T_2;
						l->node[x][y][z][12] -= T_2;
						l->node[x][y][z][14] += T_2;
					}

		pthread_barrier_wait(&barrier);

	/* Propagate */
		for(x = chunks[id]; x < chunks[id + 1]; ++x)
			for(y = 0; y < l->ly; y++)
				for(z = 0; z < l->lz; z++){
				
					//compute upper and right next neighbour nodes
					x_e = (x + 1) % l->lx;
					y_u = (y + 1) % l->ly;
					z_n = (z + 1) % l->lz;
				
					//compute lower and left next neighbour nodes
					x_w = (x - 1 + l->lx) % l->lx;
					y_d = (y - 1 + l->ly) % l->ly;
					z_s = (z - 1 + l->lz) % l->lz;
					//density propagation
					
					//zero
					l->temp[x  ][y  ][z  ][0]  = l->node[x][y][z][0];

					//east
					l->temp[x_e][y  ][z  ][1]  = l->node[x][y][z][1];
					//north
					l->temp[x_e][y_u][z  ][2]  = l->node[x][y][z][2];
					//west
					l->temp[x  ][y_u][z  ][3]  = l->node[x][y][z][3];
					//south
					l->temp[x_w][y_u][z  ][4]  = l->node[x][y][z][4];
					//down
					l->temp[x_w][y  ][z  ][5]  = l->node[x][y][z][5];
					//up
					l->temp[x_w][y_d][z  ][6]  = l->node[x][y][z][6];
				
					//east down
					l->temp[x  ][y_d][z  ][7]  = l->node[x][y][z][7];
					//east up
					l->temp[x_e][y_d][z  ][8]  = l->node[x][y][z][8];
					//north-east
					l->temp[x_e][y  ][z_n][9]  = l->node[x][y][z][9];
					//north-down
					l->temp[x  ][y  ][z_n][10] = l->node[x][y][z][10];
					//north-up
					l->temp[x_w][y  ][z_n][11] = l->node[x][y][z][11];
					//north-west
					l->temp[x_w][y  ][z_s][12] = l->node[x][y][z][12];
					//west down
					l->temp[x  ][y  ][z_s][13] = l->node[x][y][z][13];
					//west up
					l->temp[x_e][y  ][z_s][14] = l->node[x][y][z][14];
					//south-west
					l->temp[x  ][y_u][z_n][15] = l->node[x][y][z][15];
					//south down
					l->temp[x  ][y_d][z_n][16] = l->node[x][y][z][16];
					//south up
					l->temp[x  ][y_d][z_s][17] = l->node[x][y][z][17];
					//south-east
					l->temp[x  ][y_u][z_s][18] = l->node[x][y][z][18];
				}
		pthread_barrier_wait(&barrier);

	/* Bounceback */
		for(x = chunks[id]; x < chunks[id + 1]; ++x)
			for(y = 0; y < l->ly; y++)
				for(z = 0; z < l->lz; z++)
					if(l->obst[x][y][z] == true){
						//east
						l->node[x][y][z][1]  = l->temp[x][y][z][5];
						//north
						l->node[x][y][z][2]  = l->temp[x][y][z][6];
						//west
						l->node[x][y][z][3]  = l->temp[x][y][z][7];
						//south
						l->node[x][y][z][4]  = l->temp[x][y][z][8];
						//down
						l->node[x][y][z][5]  = l->temp[x][y][z][1];
						//up
						l->node[x][y][z][6]  = l->temp[x][y][z][2];
						
						//east down
						l->node[x][y][z][7]  = l->temp[x][y][z][3];
						//east up
						l->node[x][y][z][8]  = l->temp[x][y][z][4];
						//north-east
						l->node[x][y][z][9]  = l->temp[x][y][z][12];
						//north-down
						l->node[x][y][z][10] = l->temp[x][y][z][13];
						//north-up
						l->node[x][y][z][11] = l->temp[x][y][z][14];
						//north-west
						l->node[x][y][z][12] = l->temp[x][y][z][9];
						//west down
						l->node[x][y][z][13] = l->temp[x][y][z][10];
						//west up
						l->node[x][y][z][14] = l->temp[x][y][z][11];
						//south-west
						l->node[x][y][z][15] = l->temp[x][y][z][17];
						//south down
						l->node[x][y][z][16] = l->temp[x][y][z][18];
						//south up
						l->node[x][y][z][17] = l->temp[x][y][z][15];
						//south-east
						l->node[x][y][z][18] = l->temp[x][y][z][16];
					}
		pthread_barrier_wait(&barrier);

	/* Relaxation */
		for(x = chunks[id]; x < chunks[id + 1]; ++x)
			for(y = 0; y < l->ly; y++)
				for(z = 0; z < l->lz; z++)
					if(l->obst[x][y][z] == false){

						d_loc = l->temp[x][y][z][0]  + l->temp[x][y][z][1]  + l->temp[x][y][z][2]  + l->temp[x][y][z][3]  + l->temp[x][y][z][4]  + l->temp[x][y][z][5]
							  + l->temp[x][y][z][6]  + l->temp[x][y][z][7]  + l->temp[x][y][z][8]  + l->temp[x][y][z][9]  + l->temp[x][y][z][10] + l->temp[x][y][z][11]
							  + l->temp[x][y][z][12] + l->temp[x][y][z][13] + l->temp[x][y][z][14] + l->temp[x][y][z][15] + l->temp[x][y][z][16] + l->temp[x][y][z][17] 
							  + l->temp[x][y][z][18];

						//x-, y- and z- velocity components
						u_x = (l->temp[x][y][z][1] + l->temp[x][y][z][2] + l->temp[x][y][z][8] + l->temp[x][y][z][9]  + l->temp[x][y][z][14]
						      -l->temp[x][y][z][4] - l->temp[x][y][z][5] - l->temp[x][y][z][6] - l->temp[x][y][z][11] - l->temp[x][y][z][13]
							  ) / d_loc;

						u_y = (l->temp[x][y][z][2] + l->temp[x][y][z][3] + l->temp[x][y][z][4] + l->temp[x][y][z][15] + l->temp[x][y][z][18]
							  -l->temp[x][y][z][5] - l->temp[x][y][z][7] - l->temp[x][y][z][8] - l->temp[x][y][z][16] - l->temp[x][y][z][17]
							  ) / d_loc;
						
						u_z = (l->temp[x][y][z][9]  + l->temp[x][y][z][10] + l->temp[x][y][z][11] + l->temp[x][y][z][15] + l->temp[x][y][z][16]
							  -l->temp[x][y][z][12] - l->temp[x][y][z][13] - l->temp[x][y][z][14] - l->temp[x][y][z][17] - l->temp[x][y][z][18]
							  ) / d_loc;
		
						//square velocity
						u_squ = u_x * u_x + u_y * u_y + u_z * u_z;

						//n- velocity components
						//only 3 speeds would be necessary
						u_n[1]  =  u_x;
						u_n[2]  =  u_x + u_y;
						u_n[3]  =  u_y;
						u_n[4]  = -u_x + u_y;
						u_n[5]  = -u_x;
						u_n[6]  = -u_x - u_y;
						u_n[7]  = -u_y;
						u_n[8]  =  u_x - u_y;

						u_n[9]  =  u_x + u_z;
						u_n[10] =  u_z;
						u_n[11] = -u_x + u_z;
						u_n[12] = -u_x - u_z;
						u_n[13] = -u_z;
						u_n[14] =  u_x - u_z;
						u_n[15] =  u_y + u_z;
						u_n[16] = -u_y + u_z;
						u_n[17] = -u_y - u_z;
						u_n[18] =  u_y - u_z;
								  
						//zero velocity density
						n_equ[0]  = T0 * d_loc * (1.0 - u_squ   / (2.0 * CS2));

						//axis speeds: factor: t1
						n_equ[1]  = T1 * d_loc * (1.0 + u_n[1]  / CS2 + u_n[1]  * u_n[1]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[2]  = T2 * d_loc * (1.0 + u_n[2]  / CS2 + u_n[2]  * u_n[2]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[3]  = T1 * d_loc * (1.0 + u_n[3]  / CS2 + u_n[3]  * u_n[3]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[4]  = T2 * d_loc * (1.0 + u_n[4]  / CS2 + u_n[4]  * u_n[4]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[5]  = T1 * d_loc * (1.0 + u_n[5]  / CS2 + u_n[5]  * u_n[5]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[6]  = T2 * d_loc * (1.0 + u_n[6]  / CS2 + u_n[6]  * u_n[6]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						
						//diagonal speeds: factor t2
						n_equ[7]  = T1 * d_loc * (1.0 + u_n[7]  / CS2 + u_n[7]  * u_n[7]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[8]  = T2 * d_loc * (1.0 + u_n[8]  / CS2 + u_n[8]  * u_n[8]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[9]  = T2 * d_loc * (1.0 + u_n[9]  / CS2 + u_n[9]  * u_n[9]  / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[10] = T1 * d_loc * (1.0 + u_n[10] / CS2 + u_n[10] * u_n[10] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[11] = T2 * d_loc * (1.0 + u_n[11] / CS2 + u_n[11] * u_n[11] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[12] = T2 * d_loc * (1.0 + u_n[12] / CS2 + u_n[12] * u_n[12] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[13] = T1 * d_loc * (1.0 + u_n[13] / CS2 + u_n[13] * u_n[13] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[14] = T2 * d_loc * (1.0 + u_n[14] / CS2 + u_n[14] * u_n[14] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[15] = T2 * d_loc * (1.0 + u_n[15] / CS2 + u_n[15] * u_n[15] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[16] = T2 * d_loc * (1.0 + u_n[16] / CS2 + u_n[16] * u_n[16] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[17] = T2 * d_loc * (1.0 + u_n[17] / CS2 + u_n[17] * u_n[17] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
						n_equ[18] = T2 * d_loc * (1.0 + u_n[18] / CS2 + u_n[18] * u_n[18] / (2.0 * CS2 * CS2) - u_squ / (2.0 * CS2));
					
						// relaxation step
						l->node[x][y][z][0]  = l->temp[x][y][z][0]  + OMEGA * (n_equ[0]  - l->temp[x][y][z][0] );
						l->node[x][y][z][1]  = l->temp[x][y][z][1]  + OMEGA * (n_equ[1]  - l->temp[x][y][z][1] );
						l->node[x][y][z][2]  = l->temp[x][y][z][2]  + OMEGA * (n_equ[2]  - l->temp[x][y][z][2] );
						l->node[x][y][z][3]  = l->temp[x][y][z][3]  + OMEGA * (n_equ[3]  - l->temp[x][y][z][3] );
						l->node[x][y][z][4]  = l->temp[x][y][z][4]  + OMEGA * (n_equ[4]  - l->temp[x][y][z][4] );
						l->node[x][y][z][5]  = l->temp[x][y][z][5]  + OMEGA * (n_equ[5]  - l->temp[x][y][z][5] );
						l->node[x][y][z][6]  = l->temp[x][y][z][6]  + OMEGA * (n_equ[6]  - l->temp[x][y][z][6] );
						l->node[x][y][z][7]  = l->temp[x][y][z][7]  + OMEGA * (n_equ[7]  - l->temp[x][y][z][7] );
						l->node[x][y][z][8]  = l->temp[x][y][z][8]  + OMEGA * (n_equ[8]  - l->temp[x][y][z][8] );
						l->node[x][y][z][9]  = l->temp[x][y][z][9]  + OMEGA * (n_equ[9]  - l->temp[x][y][z][9] );
						l->node[x][y][z][10] = l->temp[x][y][z][10] + OMEGA * (n_equ[10] - l->temp[x][y][z][10]);
						l->node[x][y][z][11] = l->temp[x][y][z][11] + OMEGA * (n_equ[11] - l->temp[x][y][z][11]);
						l->node[x][y][z][12] = l->temp[x][y][z][12] + OMEGA * (n_equ[12] - l->temp[x][y][z][12]);
						l->node[x][y][z][13] = l->temp[x][y][z][13] + OMEGA * (n_equ[13] - l->temp[x][y][z][13]);
						l->node[x][y][z][14] = l->temp[x][y][z][14] + OMEGA * (n_equ[14] - l->temp[x][y][z][14]);
						l->node[x][y][z][15] = l->temp[x][y][z][15] + OMEGA * (n_equ[15] - l->temp[x][y][z][15]);
						l->node[x][y][z][16] = l->temp[x][y][z][16] + OMEGA * (n_equ[16] - l->temp[x][y][z][16]);
						l->node[x][y][z][17] = l->temp[x][y][z][17] + OMEGA * (n_equ[17] - l->temp[x][y][z][17]);
						l->node[x][y][z][18] = l->temp[x][y][z][18] + OMEGA * (n_equ[18] - l->temp[x][y][z][18]);
						
					}
		pthread_barrier_wait(&barrier);

	}
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
			printf("Lattice Boltzmann Method - D3Q19 - Matrix - Pthreads - Threads: %ld - Double Precision\n", num_threads);
		#else
			printf("Lattice Boltzmann Method - D3Q19 - Matrix - Pthreads - Threads: %ld - Single Precision\n", num_threads);
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

	//Threads
	long int i;
	pthread_t *thread;
	pthread_attr_t attr;
	int error;
	void *status;

	//Execution time
	double timer = {0.0};

	//Check error
	double e_i = 0.0, e_j = 0.0; /* Iteration i and j error */
	double e_r = 0.0; /* Relative error */
	double const e_s = 0.5E-10; /* Estimated error */

	/* Properties structure */
	print_parameters(argv, argc);

	/* Lattice structure */
	lattice =(s_lattice*) read_obstacles(argv, argc);

	
	init_density(lattice, checkFlag(argv, argc, "debug"));

	if(checkFlag(argv, argc, "check")){
		e_i = check_density(lattice, checkFlag(argv, argc, "debug"));
	}
	
	if(checkFlag(argv, argc, "debug")){
		printf("Start main loop\n");
	}

	/* Chunk*/
	chunks = (unsigned int *) malloc((num_threads + 1) * sizeof(unsigned int));
	for(i = 0; i < num_threads; i++)
		chunks[i] = i * (lattice->lx / num_threads);

	for(i = 1; i <= lattice->lx % num_threads; i++)
		chunks[i] += i;

	for(i = (lattice->lx % num_threads) + 1; i < num_threads; i++)
		chunks[i] += lattice->lx % num_threads;
	chunks[num_threads] = lattice->lx;


	/* Posix Threads */
	thread = (pthread_t *) malloc(num_threads * sizeof(pthread_t));

	pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);	

    pthread_barrier_init(&barrier, NULL, num_threads);	

	timer = crono();

	for(i = 0; i < num_threads; ++i){
		if(checkFlag(argv, argc, "debug"))
			printf("Creating thread %ld\n", i);
		
		error = pthread_create(&thread[i], &attr, boltzmann, (void *)i);
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

	    comp_rey(lattice, ITERATIONS, checkFlag(argv, argc, "debug"));
	    printf("End of the execution\n\n");
	}

	if(checkFlag(argv, argc, "check")){
		e_j = check_density(lattice, checkFlag(argv, argc, "debug"));
		e_r =(e_j - e_i) / e_j;

		if(e_r > e_s){
			fprintf(stderr, "Relative error is %.15f > %lg\n", e_r, e_s);
			HANDLE_ERROR();
		}
	}

	if(checkFlag(argv, argc, "out"))
		write_results(flagValueText(argv, argc, "output"), lattice);
	dealloc_lattice(lattice, checkFlag(argv, argc, "debug"));

	free(thread);
	free(lattice);
	free(chunks);

	fprintf(stderr, "%.10lf\n", timer);
	pthread_exit((void *) 0);
}