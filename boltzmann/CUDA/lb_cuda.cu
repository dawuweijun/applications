/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */

#include <stdio.h>
#include <stdlib.h>
#include "helper.cuh"
#include "lb_3D.cuh"
 #include <cuda_runtime.h>

static void GPUHandleError( cudaError_t err, const char *file, const int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

#define GPU_HANDLE_ERROR( err ) (GPUHandleError( err, __FILE__, __LINE__ ))

/* It is interesting to redistribute de forces to all points */
#ifdef USE_DOUBLE
	__global__ void redistribute(unsigned short int *obst, double *node, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#else
	__global__ void redistribute(unsigned short int *obst, float *node, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#endif
	{
		unsigned short int x = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short int y = blockIdx.y * blockDim.y + threadIdx.y;
		unsigned short int z = blockIdx.z * blockDim.z + threadIdx.z;
		
		if(x < lx && y < ly && z < lz && obst[x * ly * lz + y * lz + z] == false){				
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  += T_1;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2]  += T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4]  -= T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]  -= T_1;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6]  -= T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8]  += T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  += T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] -= T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] -= T_2;
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] += T_2;
		}
	}


//////////////////////////////////////////
// Propagate
//////////////////////////////////////////
#ifdef USE_DOUBLE
	__global__ void propagate(double *node, double *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#else
	__global__ void propagate(float *node, float *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#endif
	{
		unsigned short int x = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short int y = blockIdx.y * blockDim.y + threadIdx.y;
		unsigned short int z = blockIdx.z * blockDim.z + threadIdx.z;
				
		if(x < lx && y < ly && z < lz){
			unsigned short int x_e, x_w, y_u, y_d, z_n, z_s;

			//compute upper and right next neighbour nodes
			x_e = (x + 1) % lx;
			y_u = (y + 1) % ly;
			z_n = (z + 1) % lz;

			//compute lower and left next neighbour nodes
			x_w = (x - 1 + lx) % lx;
			y_d = (y - 1 + ly) % ly;
			z_s = (z - 1 + lz) % lz;
			//density propagation
			
			//zero
			temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0];

			//east
			temp[x_e * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1];
			//north
			temp[x_e * ly * lz * NDIM + y_u * lz * NDIM + z * NDIM + 2]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2];
			//west
			temp[x * ly * lz * NDIM + y_u * lz * NDIM + z * NDIM + 3]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3];
			//south
			temp[x_w * ly * lz * NDIM + y_u * lz * NDIM + z * NDIM + 4]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4];
			//down
			temp[x_w * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5];
			//up
			temp[x_w * ly * lz * NDIM + y_d * lz * NDIM + z * NDIM + 6]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6];

			//east down
			temp[x * ly * lz * NDIM + y_d * lz * NDIM + z * NDIM + 7]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7];
			//east up
			temp[x_e * ly * lz * NDIM + y_d * lz * NDIM + z * NDIM + 8]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8];
			//north-east
			temp[x_e * ly * lz * NDIM + y * lz * NDIM + z_n * NDIM + 9]  = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9];
			//north-down
			temp[x * ly * lz * NDIM + y * lz * NDIM + z_n * NDIM + 10] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10];
			//north-up
			temp[x_w * ly * lz * NDIM + y * lz * NDIM + z_n * NDIM + 11] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11];
			//north-west
			temp[x_w * ly * lz * NDIM + y * lz * NDIM + z_s * NDIM + 12] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12];
			//west down
			temp[x * ly * lz * NDIM + y * lz * NDIM + z_s * NDIM + 13] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13];
			//west up
			temp[x_e * ly * lz * NDIM + y * lz * NDIM + z_s * NDIM + 14] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14];
			//south-west
			temp[x * ly * lz * NDIM + y_u * lz * NDIM + z_n * NDIM + 15] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15];
			//south down
			temp[x * ly * lz * NDIM + y_d * lz * NDIM + z_n * NDIM + 16] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16];
			//south up
			temp[x * ly * lz * NDIM + y_d * lz * NDIM + z_s * NDIM + 17] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17];
			//south-east
			temp[x * ly * lz * NDIM + y_u * lz * NDIM + z_s * NDIM + 18] = node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18];
		}
	}

//////////////////////////////////////////
// Bounceback
//////////////////////////////////////////
#ifdef USE_DOUBLE
	__global__ void bounceback(unsigned short int *obst, double *node, double *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#else
	__global__ void bounceback(unsigned short int *obst, float *node, double *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#endif
	{
		unsigned short int x = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short int y = blockIdx.y * blockDim.y + threadIdx.y;
		unsigned short int z = blockIdx.z * blockDim.z + threadIdx.z;
				
		if(x < lx && y < ly && z < lz && obst[x * ly * lz + y * lz + z] == true){
			//east
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5];
			//north
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6];
			//west
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7];
			//south
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8];
			//down
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1];
			//up
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2];
			
			//east down
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3];
			//east up
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4];
			//north-east
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12];
			//north-down
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13];
			//north-up
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14];
			//north-west
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9];
			//west down
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10];
			//west up
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11];
			//south-west
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17];
			//south down
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18];
			//south up
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15];
			//south-east
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16];
		}
	}


//////////////////////////////////////////
// Relaxation
//////////////////////////////////////////
#ifdef USE_DOUBLE
	__global__ void relaxation(unsigned short int *obst, double *node, double *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#else
	__global__ void relaxation(unsigned short int *obst, float *node, double *temp, unsigned short int lx, unsigned short int ly, unsigned short int lz)
#endif
	{
		unsigned short int x = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short int y = blockIdx.y * blockDim.y + threadIdx.y;
		unsigned short int z = blockIdx.z * blockDim.z + threadIdx.z;
				
		if(x < lx && y < ly && z < lz && obst[x * ly * lz + y * lz + z] == false){
			#ifdef USE_DOUBLE
				double u_x, u_y, u_z;
				double u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
			#else
				float u_x, u_y, u_z;
				float u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
			#endif

			d_loc = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]
				  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11]
				  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17] 
				  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18];

			//x-, y- and z- velocity components
			u_x = (temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14]
			      -temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13]
				  ) / d_loc;

			u_y = (temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18]
				  -temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17]
				  ) / d_loc;
			
			u_z = (temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] + temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16]
				  -temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18]
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
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0]  + OMEGA * (n_equ[0]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 0] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1]  + OMEGA * (n_equ[1]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 1] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2]  + OMEGA * (n_equ[2]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 2] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3]  + OMEGA * (n_equ[3]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 3] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4]  + OMEGA * (n_equ[4]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 4] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5]  + OMEGA * (n_equ[5]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 5] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6]  + OMEGA * (n_equ[6]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 6] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7]  + OMEGA * (n_equ[7]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 7] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8]  + OMEGA * (n_equ[8]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 8] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9]  + OMEGA * (n_equ[9]  - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 9] );
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10] + OMEGA * (n_equ[10] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 10]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11] + OMEGA * (n_equ[11] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 11]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12] + OMEGA * (n_equ[12] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 12]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13] + OMEGA * (n_equ[13] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 13]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14] + OMEGA * (n_equ[14] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 14]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15] + OMEGA * (n_equ[15] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 15]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16] + OMEGA * (n_equ[16] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 16]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17] + OMEGA * (n_equ[17] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 17]);
			node[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18] = temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18] + OMEGA * (n_equ[18] - temp[x * ly * lz * NDIM + y * lz * NDIM + z * NDIM + 18]);	
		}
	}

int main(int argc, char **argv){

	int block;
	if(flagValueInt(argv, argc, "block") > 0){
		block = flagValueInt(argv, argc, "block");
	}
	else{
		block = 8;
	}	
	
	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Lattice Boltzmann Method - D3Q19 - Vector - CUDA - Block: %d - Double Precision\n", block);
		#else
			printf("Lattice Boltzmann Method - D3Q19 - Vector - CUDA - Block: %d - Single Precision\n", block);
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

	/* Iteration counter */
	unsigned short int time;

	//Execution time
	float timer[6] = {0.0};

	//Check error
	double e_i = 0.0, e_j = 0.0; /* Iteration i and j error */
	double e_r = 0.0; /* Relative error */
	double const e_s = 0.5E-10; /* Estimated error */

	/* Properties structure */
	print_parameters(argv, argc);

	/* Lattice structure */
	s_lattice *lattice =(s_lattice*) read_obstacles(argv, argc);

	
	init_density(lattice, checkFlag(argv, argc, "debug"));

	if(checkFlag(argv, argc, "check")){
		e_i = check_density(lattice, checkFlag(argv, argc, "debug"));
	}
	
	if(checkFlag(argv, argc, "debug")){
		printf("Start main loop\n");
	}

	// Workstation UNIPAMPA: 0 to Tesla C2075, 1 to Quadro 5000
	//GPU_HANDLE_ERROR(cudaSetDevice(0));

	// Device Memory
	unsigned short int *d_obst;
	#ifdef USE_DOUBLE
		double *d_node, *d_temp;
	#else
		float *d_node, *d_temp;
	#endif

	// Blocks and Grids.
	dim3 BLOCK_REDISTRIBUTE(block, block, block);
	dim3 BLOCK_PROPAGATE(block, block, block);
	dim3 BLOCK_BOUNCEBACK(block, block, block);
	dim3 BLOCK_RELAXATION(block, block, block);

	int a = (lattice->lx + BLOCK_REDISTRIBUTE.x - 1) / BLOCK_REDISTRIBUTE.x;
	int b = (lattice->ly + BLOCK_REDISTRIBUTE.y - 1) / BLOCK_REDISTRIBUTE.y;
	int c = (lattice->lz + BLOCK_REDISTRIBUTE.z - 1) / BLOCK_REDISTRIBUTE.z;
	dim3 GRID_REDISTRIBUTE(a, b, c);

	a = (lattice->lx + BLOCK_PROPAGATE.x - 1) / BLOCK_PROPAGATE.x;
	b = (lattice->ly + BLOCK_PROPAGATE.y - 1) / BLOCK_PROPAGATE.y;
	c = (lattice->lz + BLOCK_PROPAGATE.z - 1) / BLOCK_PROPAGATE.z;
	dim3 GRID_PROPAGATE(a, b, c);

	a = (lattice->lx + BLOCK_BOUNCEBACK.x - 1) / BLOCK_BOUNCEBACK.x;
	b = (lattice->ly + BLOCK_BOUNCEBACK.y - 1) / BLOCK_BOUNCEBACK.y;
	c = (lattice->lz + BLOCK_BOUNCEBACK.z - 1) / BLOCK_BOUNCEBACK.z;
	dim3 GRID_BOUNCEBACK(a, b, c);

	a = (lattice->lx + BLOCK_RELAXATION.x - 1) / BLOCK_RELAXATION.x;
	b = (lattice->ly + BLOCK_RELAXATION.y - 1) / BLOCK_RELAXATION.y;
	c = (lattice->lz + BLOCK_RELAXATION.z - 1) / BLOCK_RELAXATION.z;
	dim3 GRID_RELAXATION(a, b, c);
	
	// Timer
	cudaEvent_t start[2], stop[2];

	GPU_HANDLE_ERROR(cudaEventCreate(&start[0]));
	GPU_HANDLE_ERROR(cudaEventCreate(&start[1]));

	GPU_HANDLE_ERROR(cudaEventCreate(&stop[0]));
	GPU_HANDLE_ERROR(cudaEventCreate(&stop[1]));

	// Memory alloc
	GPU_HANDLE_ERROR(cudaMalloc((void **) &d_obst, lattice->lx * lattice->ly * lattice->lz * sizeof(unsigned short int)));
	#ifdef USE_DOUBLE
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(double)));
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_temp, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(double)));
	#else
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(float)));
		GPU_HANDLE_ERROR(cudaMalloc((void **) &d_temp, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(float)));
	#endif

	// Memory Copy
	GPU_HANDLE_ERROR(cudaMemcpy(d_obst, lattice->obst, lattice->lx * lattice->ly * lattice->lz * sizeof(unsigned short int), cudaMemcpyHostToDevice));
	#ifdef USE_DOUBLE
		GPU_HANDLE_ERROR(cudaMemcpy(d_node, lattice->node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(double), cudaMemcpyHostToDevice));
		GPU_HANDLE_ERROR(cudaMemcpy(d_temp, lattice->temp, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(double), cudaMemcpyHostToDevice));
	#else
		GPU_HANDLE_ERROR(cudaMemcpy(d_node, lattice->node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(float), cudaMemcpyHostToDevice));
		GPU_HANDLE_ERROR(cudaMemcpy(d_temp, lattice->temp, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(float), cudaMemcpyHostToDevice));
	#endif

	// Synchronize
	GPU_HANDLE_ERROR(cudaDeviceSynchronize());

	GPU_HANDLE_ERROR(cudaEventRecord(start[0], 0));
	for(time = 1; time <= ITERATIONS; time++){
		GPU_HANDLE_ERROR(cudaEventRecord(start[1], 0));
	    	redistribute<<<GRID_REDISTRIBUTE, BLOCK_REDISTRIBUTE>>>(d_obst, d_node, lattice->lx, lattice->ly, lattice->lz);
	    GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());
	    GPU_HANDLE_ERROR(cudaEventRecord(stop[1], 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop[1]));
  		GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer[AUX], start[1], stop[1]));
	    timer[REDISTRIBUTE] += timer[AUX];

	    GPU_HANDLE_ERROR(cudaEventRecord(start[1], 0));
	    	propagate<<<GRID_PROPAGATE, BLOCK_PROPAGATE>>>(d_node, d_temp, lattice->lx, lattice->ly, lattice->lz);
	    GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());
	    GPU_HANDLE_ERROR(cudaEventRecord(stop[1], 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop[1]));
  		GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer[AUX], start[1], stop[1]));
	    timer[PROPAGATE] += timer[AUX];

	    GPU_HANDLE_ERROR(cudaEventRecord(start[1], 0));
	    	bounceback<<<GRID_BOUNCEBACK, BLOCK_BOUNCEBACK>>>(d_obst, d_node, d_temp, lattice->lx, lattice->ly, lattice->lz);
	    GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());
	    GPU_HANDLE_ERROR(cudaEventRecord(stop[1], 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop[1]));
  		GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer[AUX], start[1], stop[1]));
	    timer[BOUNCEBACK] += timer[AUX];

	    GPU_HANDLE_ERROR(cudaEventRecord(start[1], 0));
	    	relaxation<<<GRID_RELAXATION, BLOCK_RELAXATION>>>(d_obst, d_node, d_temp, lattice->lx, lattice->ly, lattice->lz);
	    GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());
	    GPU_HANDLE_ERROR(cudaEventRecord(stop[1], 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop[1]));
  		GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer[AUX], start[1], stop[1]));
	    timer[RELAXATION] += timer[AUX];

	    if(checkFlag(argv, argc, "check")){
	    	printf("%d - ", time);
			check_density(lattice, checkFlag(argv, argc, "debug"));
		}
	}
	GPU_HANDLE_ERROR(cudaEventRecord(stop[0], 0));	GPU_HANDLE_ERROR(cudaEventSynchronize(stop[0]));
  	GPU_HANDLE_ERROR(cudaDeviceSynchronize()); GPU_HANDLE_ERROR(cudaGetLastError());
    GPU_HANDLE_ERROR(cudaEventElapsedTime(&timer[EXECUTION], start[0], stop[0]));

	// Converting miliseconds to seconds
	for(unsigned int i = 0; i < 5; ++i)
		timer[i] /= 1000;

  	// Stop events
  	GPU_HANDLE_ERROR(cudaEventDestroy(start[0]));
	GPU_HANDLE_ERROR(cudaEventDestroy(start[1]));

	GPU_HANDLE_ERROR(cudaEventDestroy(stop[0]));
	GPU_HANDLE_ERROR(cudaEventDestroy(stop[1]));

	// Copy device to host
	#ifdef USE_DOUBLE
	GPU_HANDLE_ERROR(cudaMemcpy(lattice->node, d_node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(double), cudaMemcpyDeviceToHost));
	#else
	GPU_HANDLE_ERROR(cudaMemcpy(lattice->node, d_node, lattice->lx * lattice->ly * lattice->lz * NDIM * sizeof(float), cudaMemcpyDeviceToHost));
	#endif

	// Synchronize
	GPU_HANDLE_ERROR(cudaDeviceSynchronize());

	// Free device memory
	GPU_HANDLE_ERROR(cudaFree(d_obst));
	GPU_HANDLE_ERROR(cudaFree(d_node));
	GPU_HANDLE_ERROR(cudaFree(d_temp));

	// Device Reset
	GPU_HANDLE_ERROR(cudaDeviceReset());

	if(checkFlag(argv, argc, "debug")){
	    printf("Finish main loop\n");
	    printf("Time\n");

	    printf("\tredistribute time: %.10f segs.\n", timer[REDISTRIBUTE]);

	    printf("\tpropagate time: %.10f segs.\n", timer[PROPAGATE]);

	    printf("\tbounceback time: %.10f segs.\n", timer[BOUNCEBACK]);

	    printf("\trelaxation time: %.10f segs.\n", timer[RELAXATION]);

	    printf("\texecution time: %.10f segs.\n", timer[EXECUTION]);

	    comp_rey(lattice, time-1, checkFlag(argv, argc, "debug"));
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

	free(lattice);

	fprintf(stderr, "%.10f %.10f %.10f %.10f %.10f\n", timer[REDISTRIBUTE], timer[PROPAGATE], timer[BOUNCEBACK], timer[RELAXATION], timer[EXECUTION]);
	return 0;
}
