#ifndef MKPC_SCL_DEC_H
#define MKPC_SCL_DEC_H

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <stack>
#include <array>
#include <limits>

using namespace std;

double sgn(double x);   // sign function(output sign of x)

// Min-Sum decoding function
double minsum_2(double a, double b);
double minsum_3(double a, double b, double c);

// Successive Cancellation List decoder
int assignInitialPath(stack<int> &inactivePathIndices, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel);
int clonePath(int l, vector<double> &pathMetric_LLR, stack<int> &inactivePathIndices, vector<bool> &activePath, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel);
void killPath(int l, vector<double> &pathMetric_LLR, stack<int> &inactivePathIndices, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel);
double * getArrayPointer_LLR(int lam, int l, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel);
int * getArrayPointer_C(int lam, int l, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel);
void recursivelyCalcLLR(int lam, int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel, vector<vector<int>> &P, int list_size);
void recursivelyUpdateC(int lam, int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel, vector<vector<int>> &P, int list_size);
void continuePaths_UnfrozenBit(int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<double> &pathMetric, stack<int> &inactivePathIndices, vector<int> &kernel, int list_size, vector<int *> &arrayPointer_info);
void MK_SCL_Decoder(double variance, vector<double> &y, vector<int> &I, vector<int> &x_hat, vector<int> &kernel, vector<vector<int>> &P, int list_size, vector<double> &PMR);

#endif