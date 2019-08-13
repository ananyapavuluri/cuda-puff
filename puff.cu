#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <fstream> 
#include <time.h>
#include <string>
#include <sstream>
#include <curand_kernel.h>
#include "cuda_runtime.h"

#define TMAX 100
#define length 1500000
#define trials 5
#define Num_IPR 10

using namespace std;

__global__ void puffsim(double Jleak, double IP3, double V24, double k24, double n24, double kn24, double nn24,
	double a24, double V42, double k42, double n42, double kn42, double nn42, double a42, int ah42, int Vh42, int Kh42, int BT, int k_on, int Kd, int k_off, double B_new,
	double h24_inf, double h42_inf, double m24_inf, double m42_inf, double *time, double *B, double *c, double *states, double *h42_track, double *m42_track);

__global__ void puffsim(double Jleak, double IP3, double V24, double k24, double n24, double kn24, double nn24,
	double a24, double V42, double k42, double n42, double kn42, double nn42, double a42, int ah42, int Vh42, int Kh42, int BT, int k_on, int Kd, int k_off, double B_new,
	double h24_inf, double h42_inf, double m24_inf, double m42_inf, double *time, double *B, double *c, double *states, double *h42_track, double *m42_track)
{
	printf("Launch kernel");

	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int block = blockIdx.x;

	double dt0 = 0.00010; // default time stepsize(second)

	double Jrelease = 200; // calcium released via a single IP3R channel(nM)

	double c0 = 0.1; // resting calcium concn

					 // IPR parameters

					 // state transition rates

	double q45 = 11.000;
	double q54 = 3330;
	double q12 = 1240;
	double q21 = 88;
	double q23 = 3;
	double q32 = 69;
	double q26 = 10500;
	double q62 = 4010;

	int Vs = 4000;
	int Ks = 12;
	int q24scale = 1;
	int q42scale = 1;
	int lambda_h42_scale = 1;

	int lambda_h24 = 40;
	int lambda_m24 = 100;
	int lambda_m42 = 100;

	double dt = dt0;
	double Y0_1;
	double Y0_2;
	double Jserca;
	double YK1_1;
	double YK1_2;
	double Y1_1;
	double Y1_2;
	double YK2_1;
	double YK2_2;
	double Y2_1;
	double Y2_2;
	double YK3_1;
	double YK3_2;
	double Y3_1;
	double Y3_2;
	double YK4_1;
	double YK4_2;
	double Y_new_1;
	double Y_new_2;
	double c_new;
	//double B_new;
	double YK1;
	double Y1;
	double YK2;
	double Y2;
	double YK3;
	double Y3;
	double YK4;
	double q24;
	double q42;
	double lambda_h42;

	double dt1[Num_IPR];

	double r1[Num_IPR];
	double cm[Num_IPR];
	double g[Num_IPR];
	double h24[Num_IPR];
	double h24_new[Num_IPR];
	double m24[Num_IPR];
	double m24_new[Num_IPR];
	double h42[Num_IPR];
	double h42_new[Num_IPR];
	double m42[Num_IPR];
	double m42_new[Num_IPR];


	double Ta_old[36];
	double Ta_new[36];
	for (int x = 0; x < 36; x++) {
		Ta_old[x] = 0;
		Ta_new[x] = 0;
	}

	double dtmin = dt0;
	int index = 0;

	double previous[6];
	int current[10];
	int count = 0;

	// Initial values
	c[idx * length + count] = c0;
	B[idx * length + count] = k_off * BT / (k_on*c0 + k_off);
	time[idx * length + count] = 0;
	int state[Num_IPR];
	for (int w = 0; w < Num_IPR; w++) {
		state[w] = 4;
	}
	
	curandState rndState;
	curand_init(clock64(), 13, 0, &rndState);
	for (int n = 0; n < Num_IPR; n++)
	{
		h24[n] = h24_inf;
		h24_new[n] = h24[n];
		m24[n] = m24_inf;
		m24_new[n] = m24[n];
		m42[n] = m42_inf;
		m42_new[n] = m42[n];
		h42[n] = h42_inf;
		h42_new[n] = h42[n];
		g[n] = 0;
		cm[n] = c0;
		r1[n] = curand_uniform(&rndState);
	}

	int hitpost = 0;
	
	while (time[idx * length + count] < TMAX) {
		/*int rr = idx * length + count;
		if (block == 1) {
			printf("index: %d\n", rr);
			double tm = time[idx * length + count];
			printf("time: %lf\n" ,tm);
		}*/
		
		int No = 0;
		for (int openind = 0; openind < Num_IPR; openind++) {
			if (state[openind] == 5 || state[openind] == 6) {
				No++;
			}
		}
		Y0_1 = c[idx * length + count];
		Y0_2 = B[idx * length + count];

		// Using Runge Kutta 4th order to solve ODEs for [Ca2+] and [fluo4]
		Jserca = Vs * Y0_1 / (Y0_1 + Ks);

		YK1_1 = dt * (Jrelease * No + Jleak - Jserca + (k_off*(BT - B[idx * length + count])) - (k_on* c[idx * length + count] * B[idx * length + count]));
		YK1_2 = dt * (k_off*(BT - B[idx * length + count]) - k_on * c[idx * length + count] * B[idx * length + count]);
		Y1_1 = Y0_1 + YK1_1 / 2;
		Y1_2 = Y0_2 + YK1_2 / 2;
		Jserca = Vs * Y1_1 / (Y1_1 + Ks);

		YK2_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y1_2) - k_on * Y1_1 * Y1_2);
		YK2_2 = dt * (k_off*(BT - Y1_2) - k_on * Y1_1 * Y1_2);
		Y2_1 = Y0_1 + YK2_1 / 2;
		Y2_2 = Y0_2 + YK2_2 / 2;
		Jserca = Vs * Y2_1 / (Y2_1 + Ks);

		YK3_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y2_2) - k_on * Y2_1 * Y2_2);
		YK3_2 = dt * (k_off*(BT - Y2_2) - k_on * Y2_1 * Y2_2);
		Y3_1 = Y0_1 + YK3_1;
		Y3_2 = Y0_2 + YK3_2;
		Jserca = Vs * Y3_1 / (Y3_1 + Ks);

		YK4_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y3_2) - k_on * Y3_1 * Y3_2);
		YK4_2 = dt * (k_off*(BT - Y3_2) - k_on * Y3_1 * Y3_2);

		Y_new_1 = Y0_1 + (YK1_1 + 2 * YK2_1 + 2 * YK3_1 + YK4_1) / 6;
		Y_new_2 = Y0_2 + (YK1_2 + 2 * YK2_2 + 2 * YK3_2 + YK4_2) / 6;
		c_new = Y_new_1; 		
		B_new = Y_new_2;


		for (int i = 0; i < Num_IPR; i++) {
			dt1[i] = dt;
		}

		for (int i = 0; i < Num_IPR; i++) {

			//in a closed state
			if (state[i] - 4.5 < 0) {
				cm[i] = c_new;
			}
			else if (state[i] - 4.5 > 0) {
				cm[i] = c_new + 120; // open state
			}
			else {
				cm[i] = c_new + 60;
			}

			// Updating transition rates

			m42_inf = powf(cm[i], n42) / (powf(k42, n42) + powf(cm[i], n42));
			YK1 = dt * lambda_m42*(m42_inf - m42[i]);
			Y1 = m42[i] + YK1 / 2;
			YK2 = dt * lambda_m42*(m42_inf - Y1);
			Y2 = m42[i] + YK2 / 2;
			YK3 = dt * lambda_m42*(m42_inf - Y2);
			Y3 = m42[i] + YK3;
			YK4 = dt * lambda_m42*(m42_inf - Y3);
			m42_new[i] = m42[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;

			h42_inf = powf(kn42, nn42) / (powf(kn42, nn42) + powf(cm[i], nn42));
			lambda_h42 = lambda_h42_scale * (ah42 + Vh42 * powf(cm[i], 7) / (powf(Kh42, 7) + powf(cm[i], 7)));
			YK1 = dt * lambda_h42 * (h42_inf - h42[i]);
			Y1 = h42[i] + YK1 / 2;
			YK2 = dt * lambda_h42 * (h42_inf - Y1);
			Y2 = h42[i] + YK2 / 2;
			YK3 = dt * lambda_h42 * (h42_inf - Y2);
			Y3 = h42[i] + YK3;
			YK4 = dt * lambda_h42 * (h42_inf - Y3);
			h42_new[i] = h42[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;

			m24_inf = powf(cm[i], n24) / (powf(k24, n24) + powf(cm[i], n24));
			YK1 = dt * lambda_m24*(m24_inf - m24[i]);
			Y1 = m24[i] + YK1 / 2;
			YK2 = dt * lambda_m24*(m24_inf - Y1);
			Y2 = m24[i] + YK2 / 2;
			YK3 = dt * lambda_m24*(m24_inf - Y2);
			Y3 = m24[i] + YK3;
			YK4 = dt * lambda_m24*(m24_inf - Y3);
			m24_new[i] = m24[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;

			h24_inf = powf(kn24, nn24) / (powf(kn24, nn24) + powf(cm[i], nn24));
			YK1 = dt * lambda_h24*(h24_inf - h24[i]);
			Y1 = h24[i] + YK1 / 2;
			YK2 = dt * lambda_h24*(h24_inf - Y1);
			Y2 = h24[i] + YK2 / 2;
			YK3 = dt * lambda_h24*(h24_inf - Y2);
			Y3 = h24[i] + YK3;
			YK4 = dt * lambda_h24*(h24_inf - Y3);
			h24_new[i] = h24[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;

			q24 = q24scale * (a24 + V24 * (1 - m24[i] * h24[i]));
			q42 = q42scale * (a42 + V42 * m42[i] * h42[i]);


			Ta_old[1] = q21;
			Ta_old[6] = q12;
			Ta_old[8] = q32;
			Ta_old[9] = q42;
			Ta_old[11] = q62;
			Ta_old[13] = q23;
			Ta_old[19] = q24;
			Ta_old[22] = q54;
			Ta_old[27] = q45;
			Ta_old[31] = q26;

			q24 = q24scale * (a24 + V24 * (1 - m24_new[i] * h24_new[i]));
			q42 = q42scale * (a42 + V42 * m42_new[i] * h42_new[i]);

			Ta_new[1] = q21;
			Ta_new[6] = q12;
			Ta_new[8] = q32;
			Ta_new[9] = q42;
			Ta_new[11] = q62;
			Ta_new[13] = q23;
			Ta_new[19] = q24;
			Ta_new[22] = q54;
			Ta_new[27] = q45;
			Ta_new[31] = q26;

			double g_old = g[i];
			double sum1 = 0;
			double sum2 = 0;
			int p = state[i] - 1;
			for (int x = 0; x < 6; ++x) {
				sum1 += Ta_old[x * 6 + p];
			}
			for (int y = 0; y < 6; ++y) {
				sum2 += Ta_new[y * 6 + p];
			}
			double g_new = g_old + (sum1 + sum2) / 2 * dt; 
			double epsilon = log(1 / r1[i]); //The threshold for an event to occur
			if (g_new >= epsilon) 
			{ 
				dt1[i] = (epsilon - g_old) / (g_new - g_old)*dt; // changing timestep based on threshold
			}
		}
		for (int j = 0; j < Num_IPR; ++j) {
			dt1[j] = abs(dt1[j]);
		}
		dtmin = dt1[0];
		index = 0;
		// finding minimum time step from receptors
		for (int j = 0; j < Num_IPR; ++j) {
			if (dt1[j] < dtmin) {
				dtmin = dt1[j];
				index = j;
			}
		}
		if (index == 0 && (dtmin / dt) == 1) {
			index = 100;
		}
		else {
			hitpost = 1;
		}
		dt = dtmin; // assigning a new time step value based on IP3 receptor that gave lowest time step

					// Runge Kutta 4th Order to solve ODE using the updated time step

		Jserca = Vs * Y0_1 / (Y0_1 + Ks);

		YK1_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - B[idx * length + count]) - k_on * c[idx * length + count] * B[idx * length + count]);
		YK1_2 = dt * (k_off*(BT - B[idx * length + count]) - k_on * c[idx * length + count] * B[idx * length + count]);
		Y1_1 = Y0_1 + YK1_1 / 2;
		Y1_2 = Y0_2 + YK1_2 / 2;
		Jserca = Vs * Y1_1 / (Y1_1 + Ks);

		YK2_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y1_2) - k_on * Y1_1 * Y1_2);
		YK2_2 = dt * (k_off*(BT - Y1_2) - k_on * Y1_1 * Y1_2);
		Y2_1 = Y0_1 + YK2_1 / 2;
		Y2_2 = Y0_2 + YK2_2 / 2;
		Jserca = Vs * Y2_1 / (Y2_1 + Ks);

		YK3_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y2_2) - k_on * Y2_1 * Y2_2);
		YK3_2 = dt * (k_off*(BT - Y2_2) - k_on * Y2_1 * Y2_2);
		Y3_1 = Y0_1 + YK3_1;
		Y3_2 = Y0_2 + YK3_2;
		Jserca = Vs * Y3_1 / (Y3_1 + Ks);

		YK4_1 = dt * (Jrelease*No + Jleak - Jserca + k_off * (BT - Y3_2) - k_on * Y3_1 * Y3_2);
		YK4_2 = dt * (k_off*(BT - Y3_2) - k_on * Y3_1 * Y3_2);

		Y_new_1 = Y0_1 + (YK1_1 + 2 * YK2_1 + 2 * YK3_1 + YK4_1) / 6;
		Y_new_2 = Y0_2 + (YK1_2 + 2 * YK2_2 + 2 * YK3_2 + YK4_2) / 6;
		c_new = Y_new_1;
		B_new = Y_new_2;

		//substitute for heaviside function
		for (int i = 0;i < Num_IPR;i++)
		{
			if (state[i] - 4.5 < 0)
			{
				cm[i] = c_new; 
			}
			else if (state[i] - 4.5 > 0)
			{
				cm[i] = c_new + 120; 
			}
			else if (state[i] - 4.5 == 0)
			{
				cm[i] = c_new + 60; 
			}

			// updating transition states
			m42_inf = powf(cm[i], n42) / (powf(k42, n42) + powf(cm[i], n42));
			YK1 = dt * lambda_m42*(m42_inf - m42[i]);
			Y1 = m42[i] + YK1 / 2;
			YK2 = dt * lambda_m42*(m42_inf - Y1);
			Y2 = m42[i] + YK2 / 2;
			YK3 = dt * lambda_m42*(m42_inf - Y2);
			Y3 = m42[i] + YK3;
			YK4 = dt * lambda_m42*(m42_inf - Y3);
			m42_new[i] = m42[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;
			h42_inf = powf(kn42, nn42) / (powf(kn42, nn42) + powf(cm[i], nn42));
			double tau_h42 = ah42 + Vh42 * powf(cm[i], 7) / (powf(Kh42, 7) + powf(cm[i], 7));
			YK1 = dt * tau_h42*(h42_inf - h42[i]);
			Y1 = h42[i] + YK1 / 2;
			YK2 = dt * tau_h42*(h42_inf - Y1);
			Y2 = h42[i] + YK2 / 2;
			YK3 = dt * tau_h42*(h42_inf - Y2);
			Y3 = h42[i] + YK3;
			YK4 = dt * tau_h42*(h42_inf - Y3);
			h42_new[i] = h42[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;
			m24_inf = powf(cm[i], n24) / (powf(k24, n24) + powf(cm[i], n24));
			YK1 = dt * lambda_m24*(m24_inf - m24[i]);
			Y1 = m24[i] + YK1 / 2;
			YK2 = dt * lambda_m24*(m24_inf - Y1);
			Y2 = m24[i] + YK2 / 2;
			YK3 = dt * lambda_m24*(m24_inf - Y2);
			Y3 = m24[i] + YK3;
			YK4 = dt * lambda_m24*(m24_inf - Y3);
			m24_new[i] = m24[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;
			h24_inf = powf(kn24, nn24) / (powf(kn24, nn24) + powf(cm[i], nn24));
			YK1 = dt * lambda_h24*(h24_inf - h24[i]);
			Y1 = h24[i] + YK1 / 2;
			YK2 = dt * lambda_h24*(h24_inf - Y1);
			Y2 = h24[i] + YK2 / 2;
			YK3 = dt * lambda_h24*(h24_inf - Y2);
			Y3 = h24[i] + YK3;
			YK4 = dt * lambda_h24*(h24_inf - Y3);
			h24_new[i] = h24[i] + (YK1 + 2 * YK2 + 2 * YK3 + YK4) / 6;

			q24 = q24scale * (a24 + V24 * (1 - m24[i] * h24[i]));
			q42 = q42scale * (a42 + V42 * m42[i] * h42[i]);

			Ta_old[1] = q21;
			Ta_old[6] = q12;
			Ta_old[8] = q32;
			Ta_old[9] = q42;
			Ta_old[11] = q62;
			Ta_old[13] = q23;
			Ta_old[19] = q24;
			Ta_old[22] = q54;
			Ta_old[27] = q45;
			Ta_old[31] = q26;

			q24 = q24scale * (a24 + V24 * (1 - m24_new[i] * h24_new[i]));
			q42 = q42scale * (a42 + V42 * m42_new[i] * h42_new[i]);

			Ta_new[1] = q21;
			Ta_new[6] = q12;
			Ta_new[8] = q32;
			Ta_new[9] = q42;
			Ta_new[11] = q62;
			Ta_new[13] = q23;
			Ta_new[19] = q24;
			Ta_new[22] = q54;
			Ta_new[27] = q45;
			Ta_new[31] = q26;
			int p = state[i] - 1; //p is row
			//int check = state[i];
			double sum5 = Ta_old[p];
			double sum6 = Ta_new[p];
			for (int j = 1; j < 6; j++) //j is column 0-5 i.e. all columns
			{
				sum5 = sum5 + Ta_old[j * 6 + p];
			}
			for (int j = 1; j < 6; j++)
			{
				sum6 = sum6 + Ta_new[j * 6 + p];
			}
			g[i] = g[i] + (sum5 + sum6) / 2 * dt;

		}
		
		//updating the loop
		c[idx * length + count + 1] = c_new;
		B[idx * length + count + 1] = B_new;
		//printf("B: %lf\n", B[idx * length + count]);
		time[idx * length + count + 1] = time[idx * length + count] + dt;
		count++;
		
		//updating track variables
		for (int q = 0; q < Num_IPR; q++)
		{
			h24[q] = h24_new[q];
			m24[q] = m24_new[q];
			h42[q] = h42_new[q];
			m42[q] = m42_new[q];
			h42_track[count * Num_IPR + q] = h42_new[q];
			m42_track[count  *Num_IPR + q] = m42_new[q];
		}

		if (index != 100) {
			q24 = q24scale * (a24 + V24 * (1 - m24[index] * h24[index]));
			q42 = q42scale * (a42 + V42 * m42[index] * h42[index]);

			Ta_new[1] = q12;
			Ta_new[6] = q21;
			Ta_new[8] = q23;
			Ta_new[9] = q24;
			Ta_new[11] = q26;
			Ta_new[13] = q32;
			Ta_new[19] = q42;
			Ta_new[22] = q45;
			Ta_new[27] = q54;
			Ta_new[31] = q62;

			double r2 = curand_uniform(&rndState);
			for (int k = 0; k <= 5; k++) 
			{
				int p = state[index] - 1;  
				int var = 6 * p + k;
				previous[k] = Ta_new[var]; 
			}
			int compare = 0;
			double sum1 = 0;
			double sum2 = 1;
			while ((sum1 / sum2) < r2)
			{
				sum1 = 0;
				sum2 = 0;
				for (int j = 0; j <= compare; j++)
				{
					sum1 += previous[j]; //The sum of the first 1, 2, ..., 6 elements of previous[]
				}
				for (int h = 0; h < 6; h++)
				{
					sum2 += previous[h]; //The sum of all 6 elements in previous[]
				}
				compare++;
				
				if (compare >= 6)
				{
					break;
				}
			}
			for (int y = 0; y < Num_IPR; y++)
			{
				current[y] = state[y];
			}
			current[index] = compare; 
			for (int y = 0; y < Num_IPR; y++)
			{
				//printf("state[y]: %lf\n", state[y]);
				state[y] = current[y]; 
				
			}
			r1[index] = curand_uniform(&rndState); //Assign new random number to r1 for the receptor that just changes its state
			g[index] = 0;
		}
		else {

			for (int y = 0; y < Num_IPR; y++)
			{
				state[y] = state[y]; //state does not change
			}
			dt = dt0;
		}
		
		for (int t = 0; t < Num_IPR; t++) {
			states[idx * Num_IPR * count + t] = state[t];
		}
	}

	printf("end kernel\n");
}

int main()
{

	double dt0 = 0.0001; // default time stepsize(second)

	double Jrelease = 200; // calcium released via a single IP3R channel(nM)

	double c0 = 0.1; // resting calcium concn

	//IPR parameters

		// state transition rates

	double q45 = 11.000;
	double q54 = 3330;
	double q12 = 1240;
	double q21 = 88;
	double q23 = 3;
	double q32 = 69;
	double q26 = 10500;
	double q62 = 4010;


	int Vs = 4000;
	int Ks = 12;

	double Jleak = Vs * 0.1 / (0.1 + Ks); // = 33 muM / s

	double IP3 = 0.1; // muM

	double V24 = 100;
	double k24 = 0.5490;
	double n24 = 6.3119;
	double kn24 = 96.9114;
	double nn24 = 0.0363;
	double a24 = 1 + 7.5 / (pow(IP3, 2) + 0.25);
	double V42 = 100;
	double k42 = 0.4;
	double n42 = 11.1414;
	double kn42 = 0.1703;
	double nn42 = 3.2287;
	double a42 = 1.8 * pow(IP3, 2) / (pow(IP3, 2) + 0.34);


	int ah42 = 1;
	int Vh42 = 100;
	int Kh42 = 20;

	//dye buffer
	int BT = 20; // total calmodulin buffer
	int k_on = 150; // on - rate or binding rate
	int Kd = 2; // dissociation constant (muM)
	int k_off = Kd * k_on; // dissociation rate
	double B_var = k_off * BT / (k_on * 0.1 + k_off); // resting free buffer concn
	double B_new = B_var;

	double h24_inf = pow(kn24, nn24) / (pow(kn24, nn24) + pow(c0, nn24));
	double m24_inf = pow(c0, n24) / (pow(k24, n24) + pow(c0, n24));
	double m42_inf = pow(c0, n42) / (pow(k42, n42) + pow(c0, n42));
	double h42_inf = pow(kn42, nn42) / (pow(kn42, nn42) + pow(c0, nn42));

	// sending memory to CUDA kernel

	
	

	double *host_B = new double[trials*length];
	double *B = new double[trials * length];
	cudaMalloc((void**)&B, trials *length * sizeof(double));
	cudaMemcpy(B, host_B, trials * length * sizeof(double), cudaMemcpyHostToDevice);


	double *host_c = new double[trials * length];
	double *host_time = new double[trials*length];
	double *time = new double[trials*length];

	for (int i = 0; i < (trials * length); i++) {
		host_c[i] = 0;
		host_time[i] = 0;
		host_B[i] = 0;
	}
	if (cudaSuccess != cudaMalloc((void**)&time, trials * length * sizeof(double))) {
		cout << "malloc fail" << endl;
	}

	if (cudaSuccess != cudaMemcpy(time, host_time, trials *  length * sizeof(double), cudaMemcpyHostToDevice)) {
		cout << "memcpy fail #1" << endl;
	}

	double *c = new double[trials*length];
	cudaMalloc((void**)&c, trials * length * sizeof(double));
	cudaMemcpy(c, host_c, trials*length * sizeof(double), cudaMemcpyHostToDevice);


	double *host_states = new double[trials * Num_IPR * length];
	double *states = new double[trials * Num_IPR * length];
	cudaMalloc((void**)&states,trials * Num_IPR * length * sizeof(double));
	cudaMemcpy(states, host_states, trials *Num_IPR * length * sizeof(double), cudaMemcpyHostToDevice);
	

	/*double *host_g = new double[trials * Num_IPR];
	double *g = new double[trials * Num_IPR];
	cudaMalloc((void**)&g, trials * Num_IPR * sizeof(double));
	cudaMemcpy(g, host_g, trials * Num_IPR * sizeof(double), cudaMemcpyHostToDevice); */


	/*double *host_cm = new double[Num_IPR];
	for (int i = 0; i < Num_IPR; i++) {
		host_cm[i] = 0.1;
 	}
	double *cm = new double[Num_IPR];
	cudaMalloc((void**)&cm, Num_IPR * sizeof(double));
	cudaMemcpy(cm, host_cm, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_r1 = new double[Num_IPR];
	double *r1 = new double[Num_IPR];
	cudaMalloc((void**)&r1, Num_IPR * sizeof(double));
	cudaMemcpy(r1, host_r1, Num_IPR * sizeof(double), cudaMemcpyHostToDevice); */

	double *host_h24 = new double[Num_IPR];
	for (int k = 0; k < Num_IPR; k++) {
		host_h24[k] = h24_inf;
	}
	double *h24 = new double[Num_IPR];
	cudaMalloc((void**)&h24, Num_IPR * sizeof(double));
	cudaMemcpy(h24, host_h24, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_h24_new = new double[Num_IPR];
	for (int i = 0; i < Num_IPR; i++) {
		host_h24_new[i] = host_h24[i];
	}
	double *h24_new = new double[Num_IPR];
	cudaMalloc((void**)&h24_new, Num_IPR * sizeof(double));
	cudaMemcpy(h24_new, host_h24_new, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_m24 = new double[Num_IPR];
	for (int k = 0; k < Num_IPR; k++) {
		host_m24[k] = m24_inf;
	}
	double *m24 = new double[Num_IPR];
	cudaMalloc((void**)&m24, Num_IPR * sizeof(double));
	cudaMemcpy(m24, host_m24, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_m24_new = new double[Num_IPR];
	for (int i = 0; i < Num_IPR; i++) {
		host_m24_new[i] = host_m24[i];
	}
	double *m24_new = new double[Num_IPR];
	cudaMalloc((void**)&m24_new, Num_IPR * sizeof(double));
	cudaMemcpy(m24_new, host_m24_new, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_h42 = new double[Num_IPR];
	for (int k = 0; k < Num_IPR; k++) {
		host_h42[k] = h42_inf;
	}
	double *h42 = new double[Num_IPR];
	cudaMalloc((void**)&h42, Num_IPR * sizeof(double));
	cudaMemcpy(h42, host_h42, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_h42_new = new double[Num_IPR];
	for (int i = 0; i < Num_IPR; i++) {
		host_h42_new[i] = host_h42[i];
	}
	double *h42_new = new double[Num_IPR];
	cudaMalloc((void**)&h42_new, Num_IPR * sizeof(double));
	cudaMemcpy(h42_new, host_h42_new, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_h42_track = new double[Num_IPR * length];
	double *h42_track = new double[Num_IPR * length];
	cudaMalloc((void**)&h42_track, Num_IPR * length * sizeof(double));

	double *host_m42 = new double[Num_IPR];
	for (int k = 0; k < Num_IPR; k++) {
		host_m42[k] = m42_inf;
	}
	double *m42 = new double[Num_IPR];
	cudaMalloc((void**)&m42, Num_IPR * sizeof(double));
	cudaMemcpy(m42, host_m42, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_m42_new = new double[Num_IPR];
	for (int i = 0; i < Num_IPR; i++) {
		host_m42_new[i] = host_m42[i];
	}
	double *m42_new = new double[Num_IPR];
	cudaMalloc((void**)&m42_new, Num_IPR * sizeof(double));
	cudaMemcpy(m42_new, host_m42_new, Num_IPR * sizeof(double), cudaMemcpyHostToDevice);

	double *host_m42_track = new double[Num_IPR * length];
	double *m42_track = new double[Num_IPR * length];
	cudaMalloc((void**)&m42_track, Num_IPR * length * sizeof(double));
	
	
	for (int i = 0; i < (Num_IPR); i++) {
		host_h42_track[i] = host_h42[i];
		host_m42_track[i] = host_m42[i];
	}
	cudaMemcpy(m42_track, host_m42_track, Num_IPR * length * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(h42_track, host_h42_track, Num_IPR * length * sizeof(double), cudaMemcpyHostToDevice);

	puffsim << <trials, 1 >> > (Jleak, IP3, V24, k24, n24, kn24, nn24, a24, V42, k42, n42, kn42, nn42, a42, ah42, Vh42, Kh42, BT, k_on, Kd, k_off, B_new,
		h24_inf, h42_inf, m24_inf, m42_inf, time, B, c, states,  h42_track, m42_track);
	cudaDeviceSynchronize();

	//cudaMemcpy(host_time, time, trials * length * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaSuccess != cudaMemcpy(host_time, time, trials * length * sizeof(double), cudaMemcpyDeviceToHost)) {
		cout << "memcpy fail" << endl;
	}
	if (cudaSuccess != cudaMemcpy(host_c, c, trials * length * sizeof(double), cudaMemcpyDeviceToHost)) { 
		cout << "memcpy fail calcium" << endl;
	}
	if (cudaSuccess != cudaMemcpy(host_B, B, trials * length * sizeof(double), cudaMemcpyDeviceToHost)) {
		cout << "memcpy fail B" << endl;
	}

	cudaFree(time);
	cudaFree(c);
	cudaFree(B);

	double **Calcium = new double *[length];
	for (int i = 0; i < length; i++) {
		Calcium[i] = new double[trials];
	}
	double **TIME = new double *[length];
	for (int i = 0; i < length; i++) {
		TIME[i] = new double[trials];
	}

	double **B_arr = new double *[length];
	for (int i = 0; i < length; i++) {
		B_arr[i] = new double[trials];
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < trials; j++) {
			Calcium[i][j] = 0.0000000000000;
			TIME[i][j] = 0.00000000000000;
			B_arr[i][j] = 0.00000000000000;

		}
	}


	
	int count = 0;
	for (int x = 0; x < trials; x++) {
		for (int y = 0; y < length; y++) {
			//cout << "x: " << x << endl;
			//cout << "y: " << y << endl;
			//cout << "hosttime: " << host_time[count] << endl;
			//cout << "count: " << count << endl;
			//cout << "Calcium: " << Calcium[x][y] << endl;
			
			Calcium[y][x] = host_c[count];
			TIME[y][x] = host_time[count];
			B_arr[y][x] = host_B[count];
			count++;
		}
	}

	ofstream cal;
	cal.open("Calcium.csv");
	while (cal.is_open()) {
		for (int i = 0; i < length; i++) {
			cal << Calcium[i][0];
			for (int j = 1; j < trials; j++) {
				cal << "," << Calcium[i][j];
			}
			cal << endl;
		}
		break;
	}
	cal.close();

	ofstream plottime;
	plottime.open("Time.csv");
	while (plottime.is_open()) {
		for (int i = 0; i < length; i++) {
			plottime << TIME[i][0];
			for (int j = 1; j < trials; j++) {
				plottime << "," << TIME[i][j];
			}
			plottime << endl;
		}
		break;
	}
	plottime.close();

	ofstream Bfile;
	Bfile.open("B.csv");
	while (Bfile.is_open()) {
		for (int i = 0; i < length; i++) {
			Bfile << B_arr[i][0];
			for (int j = 1; j < trials; j++) {
				Bfile << "," << B_arr[i][j];
			}
			Bfile << endl;
		}
		break;
	}
	Bfile.close();

	printf("Program finished.");
	
		//USE CUDAFREE()
	// delete pointers
	delete[] host_time;
	delete[] host_c;
	delete[] host_B;
	cudaFree(states);
	cudaFree(h42_track);
	delete[] host_h42_track;
	cudaFree(m42_track);
	delete[] host_m42_track; 

}