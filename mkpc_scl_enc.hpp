#ifndef MKPC_SCL_ENC_H
#define MKPC_SCL_ENC_H

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

using namespace std;

double Qfunc(double x); // Q function expressed by error function

void printMtx(vector<vector<int>> &M);  // output matrix

// encoding
void PassKernel(int T_size, vector<int> &u, vector<int> &x);
void Permutation(vector<int> &x, vector<int> &kernel, int lamda, vector<vector<int>> &P);
void encode(vector<int> &u, vector<int> &x, vector<int> &kernel, vector<vector<int>> &P);

// gaussian approximation
double PHI(double m);
double PHI_pr(double x);
double PHI_inv(double m);
double phi_2(double m_0, double m_1);
double phi_3(double m_0, double m_1, double m_2);
void GaussianApproximation(vector<int> &kernel, double init_value_bhattacharyya, vector<double> &bhat);

// kronecker product & constructing T_N
void vecKron_double(vector<double> &A, vector<int> &B, vector<double> &C);
void Kron(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C);
void constructT_N(vector<int> &kernel, vector<vector<int>> &T_N);

void listPartition(int k_len, int p_1, int p_2, deque<deque<int>> &K);

void recursive_comb(int *indexes, int s, int rest, function<void(int *)> f);
void foreach_comb(int n, int k, function<void(int *)> f);

// Sphere Constraint Enumeration Method
int SCEM(vector<vector<int>> &G, list<int> &R, int d_min);

// calculate minimum hamming distance from generator matrix G_k
void MtxRowSwap(vector<vector<int>> &A, int a, int b);
void MtxColSwap(vector<vector<int>> &A, int a, int b);
void gaussianElimination(vector<vector<int>> &A, int icr);
void update_w_max_and_K(int d_ub, int k, int m, vector<vector<vector<int>>> &G_sys, int w_max, vector<vector<vector<int>>> &K);
int MinDist(vector<vector<int>> &T_p);

// hybrid method for determine the position of information bits
void Spectrum_ofKroneckerProduct_ofKernels(vector<vector<int>> &R_p1, vector<vector<int>> &R_p2, vector<vector<int>> &R_p, vector<int> &S_Tp, vector<vector<int>> &T_p);
void DecideInfoSet_MinDist(list<int> &I, vector<double> &S_N, vector<vector<int>> &R_p, int K);
void hybridDesignFrozen(int K, vector<int> &kernel, int psi, double variance, vector<double> &d, vector<int> &I_0);

// calculate bhattacharyya parameter
void CalcBhattacharyya(vector<int> &kernel, double init_value_bhattacharyya, vector<double> &bhat);

#endif