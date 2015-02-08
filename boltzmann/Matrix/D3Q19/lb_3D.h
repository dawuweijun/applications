/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#ifndef LB_3D
#define LB_3D

/* Macroscopic information */
#define ITERATIONS   50    /* Maximum number of iterations */
#define DENSITY      0.1   /* Fluid density per link in g/cm³ */
#define ACCELERATION 0.05  /* Macroscopic accelleration in cm²/s */
#define OMEGA        1.85  /* Relaxation parameter */
#define REYNOLDS     10    /* Linear dimension for Reynolds number in cm */

 /* The sound velocity squared */
#define     CS2 	    (1.0/3.0)

/* The eq. coeff. has been copied from Y.H. Qian et al. Europhys. Lett. 17 (6) 479 */
#define  	T0		 	(1.0/3.0)
#define  	T1  		(1.0/18.0)
#define		T2		 	(1.0/36.0)

#define 	T_1			DENSITY * ACCELERATION * T1
#define 	T_2			DENSITY * ACCELERATION * T2

#define NDIM 19 /* D3Q19 - lattice dimension elements */


/* struct contends lattice structure*/
typedef struct {
	unsigned short int lx; /* nodes number in axis x */
	unsigned short int ly; /* nodes number in axis y */
	unsigned short int lz; /* nodes number in axis z */
	unsigned short int n_obst; /* nodes number is obstacle */
	unsigned short int ***obst; /* Obstacle Array lx * ly * lz */
	#ifdef USE_DOUBLE
		double ****node; /* n-speed lattice  n * lx * ly * lz */
		double ****temp; /* temporarily storage of fluid densities */
	#else
		float ****node; /* n-speed lattice  n * lx * ly * lz */
		float ****temp; /* temporarily storage of fluid densities */
	#endif
} s_lattice;

/* Alloc memory space to the grid */
void alloc_lattice(s_lattice *l, unsigned short int debug);

/* Returns velocity of middle element in axis x*/
#ifdef USE_DOUBLE
	double calc_velocity(s_lattice *l, unsigned short int debug);
#else
	float calc_velocity(s_lattice *l, unsigned short int debug);
#endif

/* Returns sum of densities in all nodes */
double check_density(s_lattice *l, unsigned short int debug);

/* Returns reynolds number of simulation in specific iteration(time) */
#ifdef USE_DOUBLE
	double comp_rey(s_lattice *l, int time, unsigned short int debug); 
#else
	float comp_rey(s_lattice *l, int time, unsigned short int debug); 
#endif

/* Free memory space of the grid */
void dealloc_lattice(s_lattice *l, unsigned short int debug);

/* Initializes lattice node with three density levels */
void init_density(s_lattice * l, unsigned short int debug);

/* Print parameters from properties structure */
void print_parameters(char **argv, int argc);

/* Read obstacles from file and save in lattice structure */
s_lattice* read_obstacles(char **argv, int argc);

/* Write velocities in x and y, and press in a file. If matlab is installed, save eps figures of simulation. */
void write_results(char *file, s_lattice *l);

#endif