/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */

#include <stdio.h>
#include <stdlib.h>
#include "../../../helper.h"
#include "../D3Q19/lb_3D.h"

/* It is interesting to redistribute de forces to all points */
void redistribute(s_lattice *l){
	unsigned short int x, y, z;
	
	for(x = 0; x < l->lx; x++)
		for(y = 0; y < l->ly; y++)
			for(z = 0; z < l->lz; z++)
				if(l->obst[x * l->ly * l->lz + y * l->lz + z] == false){
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  += T_1;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2]  += T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4]  -= T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]  -= T_1;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6]  -= T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8]  += T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  += T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] -= T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] -= T_2;
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] += T_2;
				}
}


//////////////////////////////////////////
// Propagate
//////////////////////////////////////////
void propagate(s_lattice *l){
	unsigned short int x, y, z;
	unsigned short int x_e, x_w, y_u, y_d, z_n, z_s;
	
	for(x = 0; x < l->lx; x++)
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
				l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0];

				//east
				l->temp[x_e * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1];
				//north
				l->temp[x_e * l->ly * l->lz * NDIM + y_u * l->lz * NDIM + z * NDIM + 2]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2];
				//west
				l->temp[x * l->ly * l->lz * NDIM + y_u * l->lz * NDIM + z * NDIM + 3]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3];
				//south
				l->temp[x_w * l->ly * l->lz * NDIM + y_u * l->lz * NDIM + z * NDIM + 4]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4];
				//down
				l->temp[x_w * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5];
				//up
				l->temp[x_w * l->ly * l->lz * NDIM + y_d * l->lz * NDIM + z * NDIM + 6]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6];
			
				//east down
				l->temp[x * l->ly * l->lz * NDIM + y_d * l->lz * NDIM + z * NDIM + 7]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7];
				//east up
				l->temp[x_e * l->ly * l->lz * NDIM + y_d * l->lz * NDIM + z * NDIM + 8]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8];
				//north-east
				l->temp[x_e * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_n * NDIM + 9]  = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9];
				//north-down
				l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_n * NDIM + 10] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10];
				//north-up
				l->temp[x_w * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_n * NDIM + 11] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11];
				//north-west
				l->temp[x_w * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_s * NDIM + 12] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12];
				//west down
				l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_s * NDIM + 13] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13];
				//west up
				l->temp[x_e * l->ly * l->lz * NDIM + y * l->lz * NDIM + z_s * NDIM + 14] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14];
				//south-west
				l->temp[x * l->ly * l->lz * NDIM + y_u * l->lz * NDIM + z_n * NDIM + 15] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15];
				//south down
				l->temp[x * l->ly * l->lz * NDIM + y_d * l->lz * NDIM + z_n * NDIM + 16] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16];
				//south up
				l->temp[x * l->ly * l->lz * NDIM + y_d * l->lz * NDIM + z_s * NDIM + 17] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17];
				//south-east
				l->temp[x * l->ly * l->lz * NDIM + y_u * l->lz * NDIM + z_s * NDIM + 18] = l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18];
			}
}

//////////////////////////////////////////
// Bounceback
//////////////////////////////////////////
void bounceback(s_lattice *l){
	unsigned short int x, y, z;

	for(x = 0; x < l->lx; x++){
		for(y = 0; y < l->ly; y++){
			for(z = 0; z < l->lz; z++){
				if(l->obst[x * l->ly * l->lz + y * l->lz + z] == true){
					//east
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5];
					//north
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6];
					//west
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7];
					//south
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8];
					//down
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1];
					//up
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2];
					
					//east down
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3];
					//east up
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4];
					//north-east
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12];
					//north-down
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13];
					//north-up
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14];
					//north-west
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9];
					//west down
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10];
					//west up
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11];
					//south-west
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17];
					//south down
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18];
					//south up
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15];
					//south-east
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16];
				}
			}
		}
	}
}


//////////////////////////////////////////
// Relaxation
//////////////////////////////////////////
void relaxation(s_lattice *l){
	//local variables
	unsigned short int x, y, z;
	#ifdef USE_DOUBLE
		double u_x, u_y, u_z;
		double u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
	#else
		float u_x, u_y, u_z;
		float u_n[NDIM], n_equ[NDIM], u_squ, d_loc;
	#endif

	for(x = 0; x < l->lx; x++){
		for(y = 0; y < l->ly; y++){
			for(z = 0; z < l->lz; z++){
				if(l->obst[x * l->ly * l->lz + y * l->lz + z] == false){

					d_loc = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]
						  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11]
						  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] 
						  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18];

					//x-, y- and z- velocity components
					u_x = (l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14]
					      -l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13]
						  ) / d_loc;

					u_y = (l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18]
						  -l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17]
						  ) / d_loc;
					
					u_z = (l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16]
						  -l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18]
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
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0]  + OMEGA * (n_equ[0]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 0] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1]  + OMEGA * (n_equ[1]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 1] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2]  + OMEGA * (n_equ[2]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 2] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3]  + OMEGA * (n_equ[3]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 3] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4]  + OMEGA * (n_equ[4]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 4] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5]  + OMEGA * (n_equ[5]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 5] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6]  + OMEGA * (n_equ[6]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 6] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7]  + OMEGA * (n_equ[7]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 7] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8]  + OMEGA * (n_equ[8]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 8] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9]  + OMEGA * (n_equ[9]  - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 9] );
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10] + OMEGA * (n_equ[10] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 10]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11] + OMEGA * (n_equ[11] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 11]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12] + OMEGA * (n_equ[12] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 12]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13] + OMEGA * (n_equ[13] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 13]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14] + OMEGA * (n_equ[14] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 14]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15] + OMEGA * (n_equ[15] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 15]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16] + OMEGA * (n_equ[16] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 16]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17] + OMEGA * (n_equ[17] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 17]);
					l->node[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18] = l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18] + OMEGA * (n_equ[18] - l->temp[x * l->ly * l->lz * NDIM + y * l->lz * NDIM + z * NDIM + 18]);
					
				}
			}
		}
	}
}

int main(int argc, char **argv){
	
	if(checkFlag(argv, argc, "debug")){
		#ifdef USE_DOUBLE
			printf("Lattice Boltzmann Method - D3Q19 - Vector -Sequential - Double Precision\n");
		#else
			printf("Lattice Boltzmann Method - D3Q19 - Vector - Sequential - Single Precision\n");
		#endif
		printf("Serpa and Schepke 2015\n");
		printf("Laboratório de Estudos Avançados - UNIPAMPA\n\n");
	}

	/* Iteration counter */
	unsigned short int time;

	//Execution time
	double timer[6] = {0.0};

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

	timer[EXECUTION] = crono();
	for(time = 1; time <= ITERATIONS; time++){
		timer[AUX] = crono();
	    redistribute(lattice);
	    timer[AUX] = crono() - timer[AUX];
	    timer[REDISTRIBUTE] += timer[AUX];

	    timer[AUX] = crono();
	    propagate(lattice);
	    timer[AUX] = crono() - timer[AUX];
	    timer[PROPAGATE] += timer[AUX];

	    timer[AUX] = crono();
	    bounceback(lattice);
	    timer[AUX] = crono() - timer[AUX];
	    timer[BOUNCEBACK] += timer[AUX];

	    timer[AUX] = crono();
	    relaxation(lattice);
	    timer[AUX] = crono() - timer[AUX];
	    timer[RELAXATION] += timer[AUX];

	    if(checkFlag(argv, argc, "check")){
	    	printf("%d - ", time);
			check_density(lattice, checkFlag(argv, argc, "debug"));
		}
	}

	timer[EXECUTION] = crono() - timer[EXECUTION];

	if(checkFlag(argv, argc, "debug")){
	    printf("Finish main loop\n");
	    printf("Time\n");

	    printf("\tredistribute time: %.10lf segs.\n", timer[REDISTRIBUTE]);

	    printf("\tpropagate time: %.10lf segs.\n", timer[PROPAGATE]);

	    printf("\tbounceback time: %.10lf segs.\n", timer[BOUNCEBACK]);

	    printf("\trelaxation time: %.10lf segs.\n", timer[RELAXATION]);

	    printf("\texecution time: %.10lf segs.\n", timer[EXECUTION]);

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

	fprintf(stderr, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", timer[REDISTRIBUTE], timer[PROPAGATE], timer[BOUNCEBACK], timer[RELAXATION], timer[EXECUTION]);
	return 0;
}
