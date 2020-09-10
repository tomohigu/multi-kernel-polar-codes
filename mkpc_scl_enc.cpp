#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <iostream>
#include <functional>
#include <list>
#include <deque>
#include <chrono>

#include "mkpc_scl_enc.hpp"

using namespace std;

#define rep(i,a,b) for(int i = a;i<b;i++)
#define REP(i,n) rep(i,0,n)
#define Q(x) erfc(x/sqrt(2.0))/2.0	// Q function

double Qfunc(double x){
	return erfc(x/sqrt(2.0))/2.0;
}

void printMtx(vector<vector<int>> &M){
	REP(i, M.size()){
		REP(j, M.at(0).size()){
			printf("%d ", M.at(i).at(j));
		}
		printf("\n");
	}
}

void PassKernel(int T_size, vector<int> &u, vector<int> &x){
	if(T_size == 2){
		REP(i, x.size()){
			if(i%2 == 0){
				x.at(i) = ( u.at(i) + u.at(i+1) ) % 2;
			}
			else{
				x.at(i) = u.at(i);
			}
		}
	}
	else if(T_size == 3){
		REP(i, x.size()){
			if(i%3 == 0){
				x.at(i) = ( u.at(i) + u.at(i+1) ) % 2;
			}
			else if(i%3 == 1){
				x.at(i) = ( u.at(i-1) + u.at(i+1) ) % 2;
			}
			else{
				x.at(i) = ( u.at(i-2) + u.at(i-1) + u.at(i) ) % 2;
			}
		}
	}
	else if(T_size == 5){
		REP(i, x.size()){
			if(i%5 == 0){
				x.at(i) = ( u.at(i) + u.at(i+1) + u.at(i+2) + u.at(i+3) ) % 2;
			}
			else if(i%5 == 1){
				x.at(i) = ( u.at(i-1) + u.at(i+2) ) % 2;
			}
			else if(i%5 == 2){
				x.at(i) = ( u.at(i-2) + u.at(i+1) + u.at(i+2) ) % 2;
			}
			else if(i%5 == 3){
				x.at(i) = ( u.at(i-3) + u.at(i-1) + u.at(i+1) ) % 2;
			}
			else{
				x.at(i) = ( u.at(i-4) + u.at(i) ) % 2;
			}
		}
	}
}

void Permutation(vector<int> &x, vector<int> &kernel, int lamda, vector<vector<int>> &P){
	vector<int> P_tmp(P.at(0).size(), 0);
	vector<int> x_tmp(x.size(), 0);
	x_tmp = x;  // copy x → x_tmp
	int N_i = x.size();
	
	if(lamda == P.size()-1){ 	// P_(m-1) = (P_0*...*P_(m-2))^(-1)
		P_tmp = P.at(0);
		rep(j, 1, P.size()-1){
			REP(i, P.at(0).size()){
				P.at(P.size()-1).at(i) = P.at(j).at(P_tmp.at(i));
			}
			P_tmp = P.at(P.size()-1);
		}
		REP(i, P.at(0).size()){
			P.at(P.size()-1).at(P_tmp.at(i)) = i;
		}
		REP(i, x.size()){
			x.at(P.at(P.size()-1).at(i)) = x_tmp.at(i);
		}
	}
	else{   // P_0, ... , P_(m-2)
		REP(i, lamda){
			N_i /= kernel.at(kernel.size()-1-i);
		}
		REP(k, x.size() / N_i){    // make Permutation
			REP(i, N_i / kernel.at(kernel.size()-1-lamda)){  // make canonical permutation Q
				REP(j, kernel.at(kernel.size()-1-lamda)){
					P.at(lamda).at(i*kernel.at(kernel.size()-1-lamda) + j + k*N_i) = j*N_i/kernel.at(kernel.size()-1-lamda) + i + k*N_i;
				}
			}
		}
		REP(i, x.size()){   // permutate x
			x.at( P.at(lamda).at(i) ) = x_tmp.at(i);
		}
	}
}

void encode(vector<int> &u, vector<int> &x, vector<int> &kernel, vector<vector<int>> &P){
	vector<int> x_tmp(x.size());
	x_tmp = u;
	REP(lamda, kernel.size()){
		PassKernel(kernel.at(kernel.size()-1-lamda), x_tmp, x);
		Permutation(x, kernel, lamda, P);
		x_tmp = x;
	}
}

double PHI(double m){
	/* if(m < 0.867861){
		return exp(0.0564*pow(m,2.0) - 0.48560*m);
	}
	else  */if(/* 0.867861 <= m && */ m < 10.0){
		return exp(-0.4527 * pow(m,0.86) + 0.0218);
	}
	else{
		return sqrt(M_PI/m) * exp(-m/4.0) * (1.0-10.0/(7.0*m));
	}
}

double PHI_pr(double x){
	if(x < 10.0){
		return -0.4527*0.86*pow(x,0.86-1.0)*exp(-0.4527*pow(x,0.86)+0.0218);
	}else{
		return -(sqrt(M_PI) * (7.0*pow(x, 2.0)+4.0*x-60.0) * exp(-x/4.0)) / (28.0*pow(x, 5.0/2.0));
	}
}

double PHI_inv(double m){
	//return pow((0.0218 - log(m))/0.4527, 1.0/0.86);
	
	// newton method
	double eps = 1.0e-15;
	double a = 0.001;
	double ah = 0.0;
	REP(i, 1000){
		ah = a - (PHI(a)-m)/PHI_pr(a);
		if(abs(ah - a) < eps){ break; }
		a = ah;
	}
	return ah;
}

double phi_2(double m_0, double m_1){
	return PHI_inv(1 - (1-PHI(m_0)) * (1-PHI(m_1)));
}

double phi_3(double m_0, double m_1, double m_2){
	return PHI_inv(1 - (1-PHI(m_0)) * (1-PHI(m_1)) * (1-PHI(m_2)));
}

void GaussianApproximation(vector<int> &kernel, double init_value_bhattacharyya, vector<double> &bhat){
	vector<vector<double>> Z(kernel.size()+1, vector<double>(bhat.size(),0));
	Z.at(0).at(0) = init_value_bhattacharyya;
	int tmp = 1;

	rep(i, 1, Z.size()){
		tmp = tmp * kernel.at(i-1);
		if(kernel.at(i-1) == 2){
			REP(j,tmp){
				if(j%2 == 0){
					Z.at(i).at(j) = phi_2(Z.at(i-1).at(j/2), Z.at(i-1).at(j/2));
				}
				else{
					Z.at(i).at(j) = 2 * Z.at(i-1).at(j/2);
				}
			}
		}
		else if(kernel.at(i-1) == 3){
			REP(j,tmp){
				if(j%3 == 0){
					Z.at(i).at(j) = phi_3(Z.at(i-1).at(j/3), Z.at(i-1).at(j/3), Z.at(i-1).at(j/3));
				}
				else if(j%3 == 1){
					Z.at(i).at(j) = Z.at(i-1).at(j/3) + phi_2(Z.at(i-1).at(j/3), Z.at(i-1).at(j/3));
				}
				else{
					Z.at(i).at(j) = 2 * Z.at(i-1).at(j/3);
				}
			}
		}
		else if(kernel.at(i-1) == 5){
			REP(j,tmp){
				if(j%5 == 0){
					Z.at(i).at(j) = phi_3(Z.at(i-1).at(j/5), Z.at(i-1).at(j/5), Z.at(i-1).at(j/5));
				}
				else if(j%5 == 1){
					Z.at(i).at(j) = phi_3(Z.at(i-1).at(j/5), Z.at(i-1).at(j/5), Z.at(i-1).at(j/5) + phi_2(Z.at(i-1).at(j/5), Z.at(i-1).at(j/5)));
				}
				else if(j%5 == 2){
					Z.at(i).at(j) = 2 * phi_2(Z.at(i-1).at(j/5), Z.at(i-1).at(j/5));
				}
				else if(j%5 == 3){
					Z.at(i).at(j) = 2 * Z.at(i-1).at(j/5) + phi_2(Z.at(i-1).at(j/5), 2 * Z.at(i-1).at(j/5));
				}
				else{
					Z.at(i).at(j) = 3 * Z.at(i-1).at(j/5);
				}
			}
		}
	}
	bhat = Z.at(kernel.size());
}

void vecKron_double(vector<double> &A, vector<int> &B, vector<double> &C){	// C[] = A ⊗ B
	REP(i, A.size()){
		REP(j, B.size()){
			C[B.size()*i + j] = A[i] * (double)B[j];
        }
    }
}

void Kron(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C){	// C[][] = A ⊗ B
	for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.at(0).size(); j++) {
			if(A.at(i).at(j) != 0){
				for (int k = 0; k < B.size(); k++) {
					for (int l = 0; l < B.at(0).size(); l++) {
						C[i*B.size() + k][j*B.at(0).size() + l] = A[i][j] * B[k][l];
					}
				}
			}
        } 
    } 
}

void constructT_N(vector<int> &kernel, vector<vector<int>> &T_N){
	vector<vector<int>> T_2 = {
		{1, 0},
		{1, 1},
	};
	vector<vector<int>> T_3 = {
		{1, 1, 1},
		{1, 0, 1},
		{0, 1, 1},
	};
	vector<vector<int>> T_5 = {
		{1, 1, 1, 1, 1},
		{1, 0, 0, 0, 0},
		{1, 0, 0, 1, 0},
		{1, 1, 1, 0, 0},
		{0, 0, 1, 1, 1},
	};
	vector<vector<int>> T_N_tmp;
	int N = 1;

	REP(i, kernel.size()){
		N *= kernel.at(i);
		if(i == 0){
			T_N_tmp.resize(N);
			REP(j, N){
				T_N_tmp.at(j).resize(N);
			}
			if(kernel.at(i) == 2){
				T_N_tmp = T_2;
			}else if(kernel.at(i) == 3){
				T_N_tmp = T_3;
			}else if(kernel.at(i) == 5){
				T_N_tmp = T_5;
			}
		}else{
			REP(j, N){
				REP(k, N){
					T_N.at(j).at(k) = 0;
				}
			}
			if(kernel.at(i) == 2){
				Kron(T_N_tmp, T_2, T_N);
			}else if(kernel.at(i) == 3){
				Kron(T_N_tmp, T_3, T_N);
			}else if(kernel.at(i) == 5){
				Kron(T_N_tmp, T_5, T_N);
			}
			T_N_tmp.resize(N);
			REP(j, N){
				T_N_tmp.at(j).resize(N);
			}
			REP(j, N){
				REP(k, N){
					T_N_tmp.at(j).at(k) = T_N.at(j).at(k);
				}
			}
		}
	}
}

void listPartition(int k_len, int p_1, int p_2, deque<deque<int>> &K){
	int i = 0;
	vector<int> intgr_pttn(k_len+1, 0);	// integer partition
	deque<int> init(p_1+1, 0);
	int k = 1;
	int l;
	int x;
	int y = k_len-1;

	while(k != 0){
		x = intgr_pttn.at(k-1) + 1;
		k -= 1;
		while(2*x <= y){
			intgr_pttn.at(k) = x;
			y -= x;
			k += 1;
		}
		l = k + 1;
		while(x <= y){
			intgr_pttn.at(k) = x;
			intgr_pttn.at(l) = y;
			intgr_pttn.at(l+1) = 0;
			if(distance(intgr_pttn.begin(), find(intgr_pttn.begin(), intgr_pttn.end(), 0)) <= p_1 && intgr_pttn.at(distance(intgr_pttn.begin(), max_element(intgr_pttn.begin(), intgr_pttn.end()))) <= p_2 ){
				if(i >= K.size()){ K.push_back(init); }
				REP(h, k+2){
					K.at(i).at(h) = intgr_pttn.at(h);
				}
				i++;
			}
			x += 1;
			y -= 1;
		}
		intgr_pttn.at(k) = x + y;
		intgr_pttn.at(k+1) = 0;
		y = x + y - 1;
		if(distance(intgr_pttn.begin(), find(intgr_pttn.begin(), intgr_pttn.end(), 0)) <= p_1 && intgr_pttn.at(distance(intgr_pttn.begin(), max_element(intgr_pttn.begin(), intgr_pttn.end()))) <= p_2 ){
			if(i >= K.size()){ K.push_back(init); }
			REP(h, k+1){
				K.at(i).at(h) = intgr_pttn.at(h);
			}
			i++;
		}
	}
}

void recursive_comb(int *indexes, int s, int rest, function<void(int *)> f) {
  if (rest == 0) {
    f(indexes);
  } else {
    if (s < 0) return;
    recursive_comb(indexes, s - 1, rest, f);
    indexes[rest - 1] = s;
    recursive_comb(indexes, s - 1, rest - 1, f);
  }
}

void foreach_comb(int n, int k, function<void(int *)> f) {
  int indexes[k];
  recursive_comb(indexes, n - 1, k, f);
}

int SCEM(vector<vector<int>> &G, list<int> &R, int d_min){
	deque<deque<int>> TAU;	// initialize
	int A_d_min = 0;	// initialize
	int i = G.size()-1;	// initialize the index of searching bit
	vector<int> v(G.at(0).size(), 0);	// initialize
	deque<int> C(G.at(0).size(), 0);	// initialize

	while(i < G.at(0).size()){
		if(any_of(R.begin(), R.end(), [&i](int x) { return x == i; })){	// information bit
			v.at(i) = 0;
			C.at(i) = 0;
			rep(j, i, v.size()){
				C.at(i) += v.at(j) * G.at(j).at(i);
			}
			C.at(i) %= 2;
			if(C.at(i) != 0){
				v.at(i) = 1;
				C.at(i) = 0;
				rep(j, i, v.size()){
					C.at(i) += v.at(j) * G.at(j).at(i);
				}
				C.at(i) %= 2;
			}
		}else{															// frozen bit
			v.at(i) = 0;
			C.at(i) = 0;
			rep(j, i, v.size()){
				C.at(i) += v.at(j) * G.at(j).at(i);
			}
			C.at(i) %= 2;
		}
	SPHERE:
		if(accumulate(C.begin()+i, C.end(), 0) <= d_min){	// satisfy the sphere constraint
			if(i > 0){
				i--;
			}else{
				TAU.push_back(C);
				A_d_min++;
				goto PRUNE;
			}
		}else{
	PRUNE:
			while(i < v.size()){
				if(any_of(R.begin(), R.end(), [&i](int x) { return x == i; }) && C.at(i) == 0){
					v.at(i) = (v.at(i) + 1) % 2;
					C.at(i) = 1;
					goto SPHERE;
				}else{
					i++;
				}
			}
		}
	}
	A_d_min--;

	/* printf("\n");
	REP(i, TAU.size()){
		REP(j, TAU.at(0).size()){
			printf("%d ", TAU.at(i).at(j));
		}
		printf("\n");
	} */

	return A_d_min;
}

void MtxRowSwap(vector<vector<int>> &A, int a, int b){	// swap row a and b
	vector<int> tmp(A.at(0).size(), 0);
	tmp = A.at(a);
	A.at(a) = A.at(b);
	A.at(b) = tmp;
}

void MtxColSwap(vector<vector<int>> &A, int a, int b){	// swap column a and b
	int tmp=0;
	REP(i, A.size()){
		tmp = A.at(i).at(a);
		A.at(i).at(a) = A.at(i).at(b);
		A.at(i).at(b) = tmp;
	}
}

void gaussianElimination(vector<vector<int>> &A, int icr){
	int lead = 0;
	int pivot;

	if(icr != 0){
		REP(i, A.size()){
			if(i+icr < A.at(0).size()){
				MtxColSwap(A, i, i+icr);
			}else{
				MtxColSwap(A, i, i+icr-A.size());
			}
		}
	}
	
	REP(r, A.size()){
		if(A.at(0).size() <= lead){ return; }
		int i = r;
		while(lead < A.at(0).size() && A.at(i).at(lead) == 0){
			i++;
			if(A.size() == i){
				i = r;
				lead++;
				if(A.at(0).size() == lead){ return; }
			}
		}
		MtxRowSwap(A, i, r);
		REP(i, A.size()){
			if(i != r){
				pivot = A.at(i).at(lead);
				REP(k, A.at(0).size()){
					A.at(i).at(k) = abs( A.at(i).at(k) - A.at(r).at(k) * pivot );
				}
			}
		}
		lead++;
	}
	
	int i=0;
	while(i<A.size() && i < A.at(0).size() && A.at(i).at(i) == 1){
		i++;
	}
	if(i != A.size()){
		int col = i;
		while(col < A.at(0).size() && A.at(i).at(col) == 0){
			col++;
		}
		MtxColSwap(A, i, col);
	}

	if(icr != 0){
		REP(i, A.size()){
			if(i+icr < A.at(0).size()){
				MtxColSwap(A, i, i+icr);
			}else{
				MtxColSwap(A, i, i+icr-A.size());
			}
		}
	}

	int rank = 0;
	while(rank<A.size() && (rank+icr < A.at(0).size() && A.at(rank).at(rank+icr) == 1 || rank+icr >= A.at(0).size() && rank+icr-A.size() < A.at(0).size() && A.at(rank).at(rank+icr-A.size()))){
		rank++;
	}
	if(rank != A.size()){
		gaussianElimination(A, icr);
	}
}

void update_w_max_and_K(int d_ub, int k, int m, vector<vector<vector<int>>> &G_sys, int w_max, vector<vector<vector<int>>> &K){
	K.assign(G_sys.begin() + m, G_sys.end());
	w_max = 1;
	int d_lb;

	do{
		d_lb = m * (w_max + 1);
		REP(j, G_sys.size()-m){
			int rank = 0;
			while(rank < G_sys.at(0).size() && rank + G_sys.at(0).size() * (m + j) < G_sys.at(0).at(0).size() && G_sys.at(m + j).at(rank).at(rank + G_sys.at(0).size() * (m + j)) == 1){
				rank++;
			}
			if(w_max >= k - rank){
				REP(i, G_sys.at(0).size()){
					REP(h, G_sys.at(0).at(0).size()){
						K.at(j).at(i).pop_back();
					}
				}
				d_lb += w_max - (k - rank) + 1;
			}
		}
		w_max++;
	}while(d_lb < d_ub && w_max <= k);
}

int MinDist(vector<vector<int>> &T_p){
	int k = T_p.size();	// the number of information bit
	int d;	// minimum distance
	vector<vector<vector<int>>> G_sys;
	int rank = 0;

	if(T_p.at(0).size() % T_p.size() == 0){		// initilization ↓
		G_sys.assign(T_p.at(0).size()/T_p.size(), T_p);
	}
	else{
		G_sys.assign(T_p.at(0).size()/T_p.size()+1, T_p);
	}
	int m = G_sys.size();
	int M = 0;
	REP(i, G_sys.size()){
		gaussianElimination(T_p, T_p.size()*i);
		G_sys.at(i) = T_p;
		rank = 0;
		while(rank < T_p.size() && rank+T_p.size()*i < T_p.at(0).size() && T_p.at(rank).at(rank+T_p.size()*i) == 1){
			rank++;
		}
		if(rank < T_p.size()){
			m = i;
			M = G_sys.size() - m;
		}
	}											// initialization ↑

	int wt_min = T_p.at(0).size();
	REP(i, G_sys.size()){
		REP(j, G_sys.at(0).size()){
			int wt = accumulate(G_sys.at(i).at(j).begin(), G_sys.at(i).at(j).end(), 0);
			if(wt < wt_min){ wt_min = wt; }
		}
	}

	int d_pr, d_ub, d_lb, w, w_max;
	d_pr = d_ub = w_max = wt_min;
	d_lb = w = 1;
	vector<vector<vector<int>>> K;
	vector<int> cw(G_sys.at(0).at(0).size(), 0);

	if(T_p.size() > 1){
		do{
			M -= K.size();
			//printf("M=%d, K_size=%d\n",M,K.size());
			//printf("comb(%d,%d)\n",k,w);
			foreach_comb(k, w, [ &k, &w, &m, &M, &G_sys, &d_pr, &d_ub, &d_lb, &w_max, &K, &cw ](int *indexes){
				REP(j, m){
					cw.assign(G_sys.at(0).at(0).size(), 0);
					REP(a, w){
						REP(b, cw.size()){
							cw.at(b) += G_sys.at(j).at(indexes[a]).at(b);
						}
					}
					REP(a, cw.size()){
						cw.at(a) %= 2;
					}
					d_pr = accumulate(cw.begin(), cw.end(), 0);
					//printf("d_pr=%d, d_ub=%d\n",d_pr,d_ub);
					if(d_pr < d_ub){
						d_ub = d_pr;
						if(d_ub <= d_lb){
							return;
						}
						update_w_max_and_K(d_ub, k, m, G_sys, w_max, K);
					}
				}
				REP(j, M){
					cw.assign(G_sys.at(0).at(0).size(), 0);
					REP(a, w){
						REP(b, cw.size()){
							cw.at(b) += G_sys.at(m+j).at(indexes[a]).at(b);
						}
					}
					REP(a, cw.size()){
						cw.at(a) %= 2;
					}
					d_pr = accumulate(cw.begin(), cw.end(), 0);
					//printf("d_pr=%d, d_ub=%d\n",d_pr,d_ub);
					if(d_pr < d_ub){
						d_ub = d_pr;
						if(d_ub <= d_lb){
							return;
						}
						update_w_max_and_K(d_ub, k, m, G_sys, w_max, K);
					}
				}
			});
			if(d_ub <= d_lb){ goto END; }

			d_lb = m * (w + 1);
			REP(j, M){
				rank = 0;
				while(rank < G_sys.at(0).size() && rank + G_sys.at(0).size() * (m + j) < G_sys.at(0).at(0).size() && G_sys.at(m + j).at(rank).at(rank + G_sys.at(0).size() * (m + j)) == 1){
					rank++;
				}
				if(w >= k - rank){
					d_lb += w - (k - rank) + 1;
				}
			}
			w++;
			//printf("d_ub=%d, d_lb=%d\n",d_ub,d_lb);
		}while(d_lb < d_ub && w <= k);
	}

	END:
	return d = d_ub;
}

void Spectrum_ofKroneckerProduct_ofKernels(vector<vector<int>> &R_p1, vector<vector<int>> &R_p2, vector<vector<int>> &R_p, vector<int> &S_Tp, vector<vector<int>> &T_p){
	vector<int> R(R_p1.at(0).size()*R_p2.at(0).size()+1, -1);

	REP(k, R_p1.size()*R_p2.size()){
		//vector<vector<int>> K(S_Tp.size(), vector<int>(R_p1.size()+1, 0));
		deque<deque<int>> K(S_Tp.size(), deque<int>(R_p1.size()+1, 0));
		listPartition(k+1, R_p1.size(), R_p2.size(), K);

		for(int l = 0; l < K.size() && K.at(l).at(0) > 0; l++){
			int idx_R = 0;
			for(int j = 0; K.at(l).at(j) > 0; j++){
				for(int i = 0; i<R_p2.at(0).size() && R_p2.at(K.at(l).at(j)-1).at(i) >= 0; i++){
					R.at(idx_R) = R_p2.at(K.at(l).at(j)-1).at(i) + R_p1.at(distance(K.at(l).begin(), find(K.at(l).begin(), K.at(l).end(), 0))-1).at(j) * R_p2.size();
					idx_R++;
				}
			}
			vector<vector<int>> G_sys(k+1, vector<int>(T_p.at(0).size()));
			for(int i=0; i < R.size() && R.at(i) >= 0; i++){
				G_sys.at(i) = T_p.at(R.at(i));
			}
			int m = MinDist(G_sys);
			if(m > S_Tp.at(k)){
				S_Tp.at(k) = m;
				REP(j, R_p1.size()*R_p2.size()){
					R_p.at(k).at(j) = R.at(j);
				}
			}
		}
	}
}

void DecideInfoSet_MinDist(list<int> &I, vector<double> &S_N, vector<vector<int>> &R_p, int K){
	REP(k, K){
		int l = distance(S_N.begin(), max_element(S_N.begin(), S_N.end()));
		int C = l % R_p.at(0).size();
		int q = (S_N.size()-l-1)/R_p.at(0).size();
		
		if(C != 0){
			for(int i=0; i < R_p.at(0).size() && R_p.at(C-1).at(i) >= 0; i++){
				auto itr = find(I.begin(), I.end(), R_p.at(C-1).at(i) + q*R_p.at(0).size());
				if(*itr == R_p.at(C-1).at(i) + q*R_p.at(0).size()){
					I.erase(itr);
				}
			}
		}
		for(int i=0; i < R_p.at(0).size() && R_p.at(C).at(i) >= 0; i++){
			I.push_front(R_p.at(C).at(i) + q * R_p.at(0).size());
			auto last_idx = I.begin();
		}
		S_N.at(l) = 0;
	}
}

void hybridDesignFrozen(int K, vector<int> &kernel, int psi, double variance, vector<double> &d, vector<int> &I_0){
	vector<double> mu;
	vector<int> S_Tp;
	vector<vector<int>> R_p;	// p = p1 * p2
	list<int> I;
	
	if(psi != 0){	// calculate reliability by GA
		vector<int> rel_kernel(psi);
		rel_kernel.assign(kernel.begin(), kernel.begin() + psi);
		int mu_size = 1;
		REP(i, rel_kernel.size()){
			mu_size *= rel_kernel.at(i);
		}
		mu.resize(mu_size);
		GaussianApproximation(rel_kernel, 2/variance, mu);
		if(psi != kernel.size()){
			reverse(mu.begin(), mu.end());
		}
	}
	if(psi != kernel.size()){	// calculate minimum distance spectrum
		vector<int> dis_kernel(kernel.size()-psi);
		dis_kernel.assign(kernel.begin() + psi, kernel.end());

		vector<vector<int>> T_2 = {
			{1, 0},
			{1, 1},
		};
		vector<int> S_T2 = {2,1};
		vector<vector<int>> R_2 = {
			{1, -1},
			{0,  1},
		};

		vector<vector<int>> T_3 = {
			{1, 1, 1},
			{1, 0, 1},
			{0, 1, 1},
		};
		vector<int> S_T3 = {3,2,1};
		vector<vector<int>> R_3 = {
			{0, -1, -1},
			{1,  2, -1},
			{0,  1,  2},
		};

		vector<vector<int>> T_5 = {
			{1, 1, 1, 1, 1},
			{1, 0, 0, 0, 0},
			{1, 0, 0, 1, 0},
			{1, 1, 1, 0, 0},
			{0, 0, 1, 1, 1},
		};
		vector<int> S_T5 = {5,3,2,1,1};
		vector<vector<int>> R_5 = {
			{0, -1, -1, -1, -1},
			{3,  4, -1, -1, -1},
			{2,  3,  4, -1, -1},
			{1,  2,  3,  4, -1},
			{0,  1,  2,  3,  4},
		};
		
		vector<vector<int>> R_p1;
		vector<vector<int>> R_p2;
		vector<vector<int>> T_p1;
		vector<vector<int>> T_p2;

		if(dis_kernel.at(0) == 2){ R_p1 = R_2;  T_p1 = T_2;  S_Tp = S_T2; }
		else if(dis_kernel.at(0) == 3){ R_p1 = R_3;  T_p1 = T_3;  S_Tp = S_T3; }
		else if(dis_kernel.at(0) == 5){ R_p1 = R_5;  T_p1 = T_5;  S_Tp = S_T5; }
		/* if(dis_kernel.at(0) == 2){ R_p1 = R_2;  S_Tp = S_T2; }
		else if(dis_kernel.at(0) == 3){ R_p1 = R_3;  S_Tp = S_T3; }
		else if(dis_kernel.at(0) == 5){ R_p1 = R_5;  S_Tp = S_T5; } */

		R_p = R_p1;	// p = p1 * p2
		
		rep(i, 1, dis_kernel.size()){
			if(dis_kernel.at(i) == 2){ R_p2 = R_2;  T_p2 = T_2; }
			else if(dis_kernel.at(i) == 3){ R_p2 = R_3;  T_p2 = T_3; }
			else if(dis_kernel.at(i) == 5){ R_p2 = R_5;  T_p2 = T_5; }
			/* if(dis_kernel.at(i) == 2){ R_p2 = R_2; }
			else if(dis_kernel.at(i) == 3){ R_p2 = R_3; }
			else if(dis_kernel.at(i) == 5){ R_p2 = R_5; } */
			vector<vector<int>> T_p(T_p1.size()*T_p2.size(), vector<int>(T_p1.at(0).size()*T_p2.at(0).size(), 0));
			Kron(T_p1, T_p2, T_p);
			S_Tp.assign(R_p1.size()*R_p2.size(), 0);
			R_p.assign(R_p1.size()*R_p2.size(), vector<int>(R_p1.at(0).size()*R_p2.at(0).size(), -1));
			Spectrum_ofKroneckerProduct_ofKernels(R_p1, R_p2, R_p, S_Tp, T_p);
			R_p1 = R_p;
			T_p1 = T_p;
		}
	}

	if(psi != 0 && psi != kernel.size()){	// hybrid design
		vecKron_double(mu, S_Tp, d);
		DecideInfoSet_MinDist(I, d, R_p, K);
	}
	else if(psi == 0){	// only distance design
		REP(i, d.size()){	// copy S_Tp → d
			d.at(i) = (double)S_Tp.at(i);
		}
		DecideInfoSet_MinDist(I, d, R_p, K);
	}
	else if(psi == kernel.size()){	// only reliability design
		I.assign(K, -1);
		for(auto itr = I.begin(); itr != I.end(); ++itr){
			*itr = distance(mu.begin(), max_element(mu.begin(), mu.end()));
			mu.at(*itr) = -1;
		}
	}

	int count = 0;
	for(auto itr = I.begin(); itr != I.end(); ++itr){
		I_0.at(count) = *itr;
		count++;
	}
}

void CalcBhattacharyya(vector<int> &kernel, double init_value_bhattacharyya, vector<double> &bhat){
	vector<vector<double>> Z(kernel.size()+1, vector<double>(bhat.size(),0));
	Z.at(0).at(0) = init_value_bhattacharyya;
	int tmp = 1;

	rep(i, 1, Z.size()){
		tmp = tmp * kernel.at(i-1);
		if(kernel.at(i-1) == 2){
			REP(j,tmp){
				if(j%2 == 0){
					Z.at(i).at(j) = 2*Z.at(i-1).at(j/2)-pow(Z.at(i-1).at(j/2),2);
				}
				else{
					Z.at(i).at(j) = pow(Z.at(i-1).at(j/2),2);
				}
			}
		}
		else if(kernel.at(i-1) == 3){
			REP(j,tmp){
				if(j%3 == 0){
					Z.at(i).at(j) = 3*Z.at(i-1).at(j/3) - 3*pow(Z.at(i-1).at(j/3),2) + pow(Z.at(i-1).at(j/3),3);
				}
				else if(j%3 == 1){
					Z.at(i).at(j) = 2*pow(Z.at(i-1).at(j/3),2) - pow(Z.at(i-1).at(j/3),3);
				}
				else{
					Z.at(i).at(j) = pow(Z.at(i-1).at(j/3),2);
				}
			}
		}
	}
	REP(i,bhat.size()){
		bhat.at(i) = -Z.at(kernel.size()).at(i);
	}
}