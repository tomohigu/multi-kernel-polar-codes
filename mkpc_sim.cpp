#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <list>
#include <deque>
#include <chrono>
#include <stack>
#include <array>
#include <limits>

#include "mkpc_scl_enc.hpp"
#include "mkpc_scl_dec.hpp"

using namespace std;

#define N_sim 1000
#define L 8	// list size

#define rep(i,a,b) for(int i = a;i<b;i++)
#define REP(i,n) rep(i,0,n)

// rand settings
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
std::uniform_real_distribution<double> rand01(0.0, 1.0); // uniform distribution [0, 1]

void make_signal(vector<int> &info){
	REP(i,info.size()) {
		if (rand01(mt) > 0.5) {
			info.at(i) = 1;
		}
		else {
			info.at(i) = 0;
		}
	}
}

void AWGN(vector<int> &x, vector<double> &y, double variance){
	std::normal_distribution<> gauss_noise(0.0, sqrt(variance) );
	REP(i,x.size()){
		if (x.at(i) == 0) {
			y.at(i) = 1.0 + gauss_noise(mt);
		}
		else {
			y.at(i) = -1.0 + gauss_noise(mt);
		}
	}
}

/* void RFChannel(vector<int> &x, vector<double> &y, double variance, vector<double> &h){
	std::normal_distribution<> gauss_noise(0.0, sqrt(variance) );
	REP(i,x.size()){
		h.at(i) = sqrt(pow(gauss_noise(mt),2) + pow(gauss_noise(mt),2));
		y.at(i) = h.at(i) * pow(-1.0,x.at(i)) + gauss_noise(mt);
	}
} */

double calcPW(int x, vector<int> &kernel){
	vector<double> beta_2 = {0, 1.1892};
	vector<double> beta_3 = {0, 1.27, 1.33 };	// only T_3
	//vector<double> beta_3 = {0, 1.8, 2.0 };	// Multi-Kernel
    double res = 0;

    // decimal 2 binal
    for (int i = 0; x>0; i++){
		if(kernel.at(i) == 2){
			res += ( (x%kernel.at(i)>0) ? 1:0 ) * pow(beta_2[x%kernel.at(i)], i);
		}
		else if(kernel.at(i) == 3){
			//res += ( (x%kernel.at(i)>0) ? 1:0 ) *   pow(beta_3[x%kernel.at(i)], i);
			if(i != 0){
				res += ( (x%kernel.at(i)>0) ? 1:0 ) * pow(beta_3[x%kernel.at(i)], i);
			}
			else{
				res += ( (x%kernel.at(i)>0) ? 1:0 ) * pow(beta_3[x%kernel.at(i)], i) + ( (x%kernel.at(i)==2) ? 0.4:0);
			}
		}
		x = x/kernel.at(i);
    }

    return res;
}

void beta_expansion(vector<int> &kernel, vector<double> &pw){
	REP(i, pw.size()){
		pw[i] =  calcPW(i, kernel);
	}
}

int main(void){
	/* FILE *fp;
	const char *fname = "HYB_729_psi_3.txt";
	fp = fopen( fname, "w" ); */

	vector<int> kernel = { 2,2,2,5,5 };
	int f = 0;
	//int e = 0;
	int N = 1;
	REP(i, kernel.size()){
		N *= kernel.at(i);
	}
	int k = N/2; // N of information bit
	double Rate = (double)k/(double)N;  // encoding rate (Rate = k/N)
	int psi = 2; // parameter determining distance and reliability design ratio ( ψ=0：distance design, ψ=kernel.size()：reliability design, otherwise：hybrid design )
	vector<int> I(k, 0);	// indeices of information bits
	

	for(double Eb_N0 = 1.0; Eb_N0 <= 5.0; Eb_N0 += 0.5){
		double variance = pow(10.0, -Eb_N0/10.0) / (2.0*Rate);
		double Z_AWGN = 2.0/variance;
		vector<double> d(N, 0.0);    // vector for sorted bhattacharyya parameter
		if(psi != 0){ I.assign(k, 0); }

		//CalcBhattacharyya(kernel, exp(-pow(10.0, Rate*Eb_N0/10.0)), d);
		//CalcBhattacharyya(kernel, exp(-1.0/(2.0*variance)), d);
		GaussianApproximation(kernel, Z_AWGN, d);
		//beta_expansion(kernel, d);
		/* if(psi == 0 && Eb_N0 == 1.0 || psi != 0){
			hybridDesignFrozen(k, kernel, psi, variance, d, I);
		} */


		// please comment out in case of hybrid design
		REP(i, k){
			int max_idx = distance(d.begin(), max_element(d.begin(), d.end()));
			I.at(i) = max_idx;
			d.at(max_idx) = -1000000;
		}
		sort(I.begin(), I.end());

		REP(i, N_sim){
			vector<int> info(k,0),u(N,0),x(N,0),u_hat(N,0);
			vector<vector<int>> P(kernel.size(), vector<int>(u.size(), 0)); // permutation matrix P
			vector<double> y(N, 0)/* , h(N, 0) */;	// channel output(y), channel condition(h)
			vector<double> PMR(N, 0.0);

			make_signal(info);                         // generate information bits
			REP(j, info.size()){ u.at(I.at(j)) = info.at(j); }
			
			encode(u, x, kernel, P); // encode x = uBG
			/* REP(j, P.size()){
				REP(m, P.at(0).size()){
					printf("%d ",P.at(j).at(m));
				}
				printf("\n");
			}
			printf("u: ");
			REP(j, u.size()){
				printf("%d ", u.at(j));
			}printf("\n");
			printf("x: ");
			REP(j, x.size()){
				printf("%d ", x.at(j));
			}printf("\n"); */
			AWGN(x,y,variance);
			MK_SCL_Decoder(variance, y, I, u_hat, kernel, P, L, PMR);

			/* printf("y: ");
			REP(j, y.size()){
				printf("%lf ", y.at(j));
			}printf("\n");
			printf("x: ");
			REP(j, u_hat.size()){
				printf("%d ", u_hat.at(j));
			}printf("\n"); */

			/* REP(j, I.size()){   // calc BER
				if(u_hat.at(I.at(j)) != x.at(I.at(j))){
					e++;
				}
			} */
			REP(j, I.size()){   // calc BLER
				if(u_hat.at(I.at(j)) != u.at(I.at(j))){
					f++;
					break;
				}
			}
		}
		//fprintf(fp, "%.2lf,%.10lf,%.10lf\n", Eb_N0, (double)f/(double)N_sim, (double)e/(double)N_sim/(double)k);
		//fprintf(fp, "%.2lf\t%.6lf\n",Eb_N0,(double)f/(double)N_sim);
		printf("%.2lf %.6lf\n", Eb_N0, (double)f/(double)N_sim);
		//e = 0;
		f = 0;
	}

	//fclose(fp);

	return 0;
}