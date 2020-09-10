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

#include "mkpc_scl_dec.hpp"

using namespace std;

#define rep(i,a,b) for(int i = a;i<b;i++)
#define REP(i,n) rep(i,0,n)
#define arctanh(x) 1.0/2.0*(log(1.0+x)-log(1.0-x))

double sgn(double x){
	return (x>0)-(x<0);
}

double minsum_2(double a, double b){
	return sgn(a) * sgn(b) * min(abs(a), abs(b));
}

double minsum_3(double a, double b, double c){
	return sgn(a) * sgn(b) * sgn(c) * min( {abs(a), abs(b), abs(c)} );
}

int assignInitialPath(stack<int> &inactivePathIndices, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel){
	int l = inactivePathIndices.top();
	inactivePathIndices.pop();
	activePath.at(l) = true;

	REP(lam, kernel.size()+1){
		int s = inactiveArrayIndices.at(lam).top();
		inactiveArrayIndices.at(lam).pop();
		pathIndexToArrayIndex.at(lam).at(l) = s;
		arrayReferenceCount.at(lam).at(s) = 1;
	}

	return l;
}

int clonePath(int l, vector<double> &pathMetric, stack<int> &inactivePathIndices, vector<bool> &activePath, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel){
	int lp = inactivePathIndices.top();
	inactivePathIndices.pop();
  	activePath.at(lp) = true;
  
	pathMetric.at(lp) = pathMetric.at(l);

	REP(lam, kernel.size()+1){
		int s = pathIndexToArrayIndex.at(lam).at(l);
		pathIndexToArrayIndex.at(lam).at(lp) = s;
		arrayReferenceCount.at(lam).at(s)++;
	}

	return lp;
}

void killPath(int l, vector<double> &pathMetric, stack<int> &inactivePathIndices, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<int> &kernel){
	activePath.at(l) = false;
	inactivePathIndices.push(l);

	pathMetric.at(l) = 0;

	REP(lam, kernel.size()+1){
		int s = pathIndexToArrayIndex.at(lam).at(l);
		arrayReferenceCount.at(lam).at(s)--;
		if(arrayReferenceCount.at(lam).at(s) == 0){
			inactiveArrayIndices.at(lam).push(s);
		}
	}
}

double * getArrayPointer_LLR(int lam, int l, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel){
	int s = pathIndexToArrayIndex.at(lam).at(l);
	int sp;
	int N_p = 1;
	REP(i, kernel.size()-lam){
		N_p *= kernel.at(kernel.size() - 1 - i);
	}

	if(arrayReferenceCount.at(lam).at(s) == 1){ sp = s; }
	else{
		sp = inactiveArrayIndices.at(lam).top();
		inactiveArrayIndices.at(lam).pop();

		// copy
		if(lam != 0 && kernel.at(lam-1) == 2){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 2*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam != 0 && kernel.at(lam-1) == 3){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 3*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam != 0 && kernel.at(lam-1) == 5){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 5*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam == 0){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + N_p , arrayPointer_C.at(lam).at(sp));
		}
		copy(arrayPointer_LLR.at(lam).at(s), arrayPointer_LLR.at(lam).at(s) + N_p , arrayPointer_LLR.at(lam).at(sp));

		arrayReferenceCount.at(lam).at(s)--;
		arrayReferenceCount.at(lam).at(sp) = 1;
		pathIndexToArrayIndex.at(lam).at(l) = sp;
	}

	return arrayPointer_LLR.at(lam).at(sp);
}

int * getArrayPointer_C(int lam, int l, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel){
	int s = pathIndexToArrayIndex.at(lam).at(l);
	int sp;
	int N_p = 1;
	REP(i, kernel.size()-lam){
		N_p *= kernel.at(kernel.size() - 1 - i);
	}

	if(arrayReferenceCount.at(lam).at(s) == 1){ sp = s; }
	else{
		sp = inactiveArrayIndices.at(lam).top();
		inactiveArrayIndices.at(lam).pop();

		// copy
		if(lam != 0 && kernel.at(lam-1) == 2){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 2*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam != 0 && kernel.at(lam-1) == 3){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 3*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam != 0 && kernel.at(lam-1) == 5){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + 5*N_p , arrayPointer_C.at(lam).at(sp));
		}else if(lam == 0){
			copy(arrayPointer_C.at(lam).at(s), arrayPointer_C.at(lam).at(s) + N_p , arrayPointer_C.at(lam).at(sp));
		}
		copy(arrayPointer_LLR.at(lam).at(s), arrayPointer_LLR.at(lam).at(s) + N_p , arrayPointer_LLR.at(lam).at(sp));

		arrayReferenceCount.at(lam).at(s)--;
		arrayReferenceCount.at(lam).at(sp) = 1;
		pathIndexToArrayIndex.at(lam).at(l) = sp;
	}

	return arrayPointer_C.at(lam).at(sp);
}

void recursivelyCalcLLR(int lam, int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel, vector<vector<int>> &P, int list_size){
	if(lam == 0){ return; }
	int N_p = 1;
	REP(i, kernel.size()-lam){
		N_p *= kernel.at(kernel.size() - 1 - i);
	}

	if(kernel.at(lam-1) == 2){
		int psi = phi/2;

		if(phi%2 == 0){ recursivelyCalcLLR(lam-1, psi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size); }

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			double * LLR_lam   = getArrayPointer_LLR(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			double * LLR_lam_1 = getArrayPointer_LLR(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int * c_lam = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					if(phi%2 == 0){
						if( 10 > max(abs(LLR_lam_1[2*Beta]), abs(LLR_lam_1[2*Beta+1])) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[2*Beta]/2.0) * tanh(LLR_lam_1[2*Beta+1]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_2(LLR_lam_1[2*Beta], LLR_lam_1[2*Beta+1]);
						}
					}
					else{
						int u0 = c_lam[2*Beta];
						LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[2*Beta] + LLR_lam_1[2*Beta+1];
					}
				}
				else{
					if(phi%2 == 0){
						if( 10 > max(abs(LLR_lam_1[P.at(P.size()-1).at(2*Beta)]), abs(LLR_lam_1[P.at(P.size()-1).at(2*Beta+1)])) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(2*Beta)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(2*Beta+1)]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_2(LLR_lam_1[P.at(P.size()-1).at(2*Beta)], LLR_lam_1[P.at(P.size()-1).at(2*Beta+1)]);
						}
					}
					else{
						int u0 = c_lam[2*Beta];
						LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[P.at(P.size()-1).at(2*Beta)] + LLR_lam_1[P.at(P.size()-1).at(2*Beta+1)];
					}
				}
				//printf("llr_%d[%d]=%lf\n", lam, Beta, LLR_lam[Beta]);
			}
		}
	}
	else if(kernel.at(lam-1) == 3){
		int psi = phi/3;

		if(phi%3 == 0){ recursivelyCalcLLR(lam-1, psi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size); }

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			double * LLR_lam   = getArrayPointer_LLR(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			double * LLR_lam_1 = getArrayPointer_LLR(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int * c_lam = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					if(phi%3 == 0){
						if( 10 > max( { abs(LLR_lam_1[3*Beta]), abs(LLR_lam_1[3*Beta+1]), abs(LLR_lam_1[3*Beta+2]) } ) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[3*Beta]/2.0) * tanh(LLR_lam_1[3*Beta+1]/2.0) * tanh(LLR_lam_1[3*Beta+2]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_3(LLR_lam_1[3*Beta], LLR_lam_1[3*Beta+1], LLR_lam_1[3*Beta+2]);
						}
					}
					else if(phi%3 == 1){
						int u0 = c_lam[3*Beta];
						if( 10 > max( abs(LLR_lam_1[3*Beta+1]), abs(LLR_lam_1[3*Beta+2]) ) ){
							LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[3*Beta] + 2.0 * arctanh( tanh(LLR_lam_1[3*Beta+1]/2.0) * tanh(LLR_lam_1[3*Beta+2]/2.0) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[3*Beta] + minsum_2(LLR_lam_1[3*Beta+1], LLR_lam_1[3*Beta+2]);
						}
					}
					else{
						int u0 = c_lam[3*Beta];
						int u1 = c_lam[3*Beta+1];
						LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[3*Beta+1] + pow(-1, (u0+u1)%2) * LLR_lam_1[3*Beta+2];
					}
				}
				else{
					if(phi%3 == 0){
						if( 10 > max( { abs(LLR_lam_1[P.at(P.size()-1).at(3*Beta)]), abs(LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)]), abs(LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]) } ) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(3*Beta)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_3(LLR_lam_1[P.at(P.size()-1).at(3*Beta)], LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)], LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]);
						}
					}
					else if(phi%3 == 1){
						int u0 = c_lam[3*Beta];
						if( 10 > max( abs(LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)]), abs(LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]) ) ){
							LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[P.at(P.size()-1).at(3*Beta)] + 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]/2.0) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[P.at(P.size()-1).at(3*Beta)] + minsum_2(LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)], LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)]);
						}
					}
					else{
						int u0 = c_lam[3*Beta];
						int u1 = c_lam[3*Beta+1];
						LLR_lam[Beta] = pow(-1, u0) * LLR_lam_1[P.at(P.size()-1).at(3*Beta+1)] + pow(-1, (u0+u1)%2) * LLR_lam_1[P.at(P.size()-1).at(3*Beta+2)];
					}
				}
				//printf("llr_%d[%d]=%lf\n", lam, Beta, LLR_lam[Beta]);
			}
		}
	}
	else if(kernel.at(lam-1) == 5){
		int psi = phi/5;

		if(phi%5 == 0){ recursivelyCalcLLR(lam-1, psi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size); }

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			double * LLR_lam   = getArrayPointer_LLR(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			double * LLR_lam_1 = getArrayPointer_LLR(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int * c_lam = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					if(phi%5 == 0){
						if( 10 > max( { abs(LLR_lam_1[5*Beta+1]), abs(LLR_lam_1[5*Beta+2]), abs(LLR_lam_1[5*Beta+4]) } ) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[5*Beta+1]/2.0) * tanh(LLR_lam_1[5*Beta+2]/2.0) * tanh(LLR_lam_1[5*Beta+4]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_3(LLR_lam_1[5*Beta+1], LLR_lam_1[5*Beta+2], LLR_lam_1[5*Beta+4]);
						}
					}
					else if(phi%5 == 1){
						int u0 = c_lam[5*Beta];
						if( 10 > max( { abs(LLR_lam_1[5*Beta]), abs(LLR_lam_1[5*Beta+1]), abs(LLR_lam_1[5*Beta+2]), abs(LLR_lam_1[5*Beta+3]), abs(LLR_lam_1[5*Beta+4]) } ) ){
							LLR_lam[Beta] = pow(-1, u0) * 2.0 * arctanh( tanh(LLR_lam_1[5*Beta]/2.0) * tanh(LLR_lam_1[5*Beta+3]/2.0) * tanh( (LLR_lam_1[5*Beta+2] + 2.0 * arctanh(tanh(LLR_lam_1[5*Beta+1]/2.0) * tanh(LLR_lam_1[5*Beta+4]/2.0))) / 2.0 ) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u0) * minsum_3( LLR_lam_1[5*Beta], LLR_lam_1[5*Beta+3], LLR_lam_1[5*Beta+2] + minsum_2(LLR_lam_1[5*Beta+1], LLR_lam_1[5*Beta+4]) );
						}
					}
					else if(phi%5 == 2){
						int u1 = c_lam[5*Beta+1];
						if( 10 > max( { abs(LLR_lam_1[5*Beta]), abs(LLR_lam_1[5*Beta+1]), abs(LLR_lam_1[5*Beta+3]), abs(LLR_lam_1[5*Beta+4]) } ) ){
							LLR_lam[Beta] = pow(-1, u1) * 2.0 * arctanh( tanh(LLR_lam_1[5*Beta]/2.0) * tanh(LLR_lam_1[5*Beta+1]/2.0) ) + 2.0 * arctanh( tanh(LLR_lam_1[5*Beta+3]/2.0) * tanh(LLR_lam_1[5*Beta+4]/2.0) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u1) * minsum_2(LLR_lam_1[5*Beta], LLR_lam_1[5*Beta+1]) + minsum_2(LLR_lam_1[5*Beta+3], LLR_lam_1[5*Beta+4]);
						}
					}
					else if(phi%5 == 3){
						int u0 = c_lam[5*Beta];
						int u1 = c_lam[5*Beta+1];
						int u2 = c_lam[5*Beta+2];
						if( 10 > max( abs(LLR_lam_1[5*Beta+2]), abs(LLR_lam_1[5*Beta+3] + LLR_lam_1[5*Beta+4]) ) ){
							LLR_lam[Beta] = pow(-1, (u0+u1+u2)%2)*LLR_lam_1[5*Beta] + pow(-1, u0)*LLR_lam_1[5*Beta+1] + 2.0 * arctanh( tanh(LLR_lam_1[5*Beta+2]/2.0) * tanh( (LLR_lam_1[5*Beta+3]+LLR_lam_1[5*Beta+4]) / 2.0 ) );
						}
						else{
							LLR_lam[Beta] = pow(-1, (u0+u1+u2)%2)*LLR_lam_1[5*Beta] + pow(-1, u0)*LLR_lam_1[5*Beta+1] + minsum_2(LLR_lam_1[5*Beta+2], LLR_lam_1[5*Beta+3] + LLR_lam_1[5*Beta+4]);
						}
					}
					else{
						int u0 = c_lam[5*Beta];
						int u2 = c_lam[5*Beta+2];
						int u3 = c_lam[5*Beta+3];
						LLR_lam[Beta] = pow(-1, (u0+u3)%2)*LLR_lam_1[5*Beta+2] + pow(-1, (u0+u2)%2)*LLR_lam_1[5*Beta+3] + pow(-1, u0)*LLR_lam_1[5*Beta+4];
					}
				}
				else{
					if(phi%5 == 0){
						if( 10 > max( { abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) } ) ){
							LLR_lam[Beta] = 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]/2.0) );
						}
						else{
							LLR_lam[Beta] = minsum_3(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]);
						}
					}
					else if(phi%5 == 1){
						int u0 = c_lam[5*Beta];
						if( 10 > max( { abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) } ) ){
							LLR_lam[Beta] = pow(-1, u0) * 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)]/2.0) * tanh( (LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)] + 2.0 * arctanh(tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]/2.0))) / 2.0 ) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u0) * minsum_3( LLR_lam_1[P.at(P.size()-1).at(5*Beta)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)] + minsum_2(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) );
						}
					}
					else if(phi%5 == 2){
						int u1 = c_lam[5*Beta+1];
						if( 10 > max( { abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) } ) ){
							LLR_lam[Beta] = pow(-1, u1) * 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]/2.0) ) + 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)]/2.0) * tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]/2.0) );
						}
						else{
							LLR_lam[Beta] = pow(-1, u1) * minsum_2(LLR_lam_1[P.at(P.size()-1).at(5*Beta)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)]) + minsum_2(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]);
						}
					}
					else if(phi%5 == 3){
						int u0 = c_lam[5*Beta];
						int u1 = c_lam[5*Beta+1];
						int u2 = c_lam[5*Beta+2];
						if( 10 > max( abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)]), abs(LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)] + LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) ) ){
							LLR_lam[Beta] = pow(-1, (u0+u1+u2)%2)*LLR_lam_1[P.at(P.size()-1).at(5*Beta)] + pow(-1, u0)*LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)] + 2.0 * arctanh( tanh(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)]/2.0) * tanh( (LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)]+LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]) / 2.0 ) );
						}
						else{
							LLR_lam[Beta] = pow(-1, (u0+u1+u2)%2)*LLR_lam_1[P.at(P.size()-1).at(5*Beta)] + pow(-1, u0)*LLR_lam_1[P.at(P.size()-1).at(5*Beta+1)] + minsum_2(LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)], LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)] + LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)]);
						}
					}
					else{
						int u0 = c_lam[5*Beta];
						int u2 = c_lam[5*Beta+2];
						int u3 = c_lam[5*Beta+3];
						LLR_lam[Beta] = pow(-1, (u0+u3)%2)*LLR_lam_1[P.at(P.size()-1).at(5*Beta+2)] + pow(-1, (u0+u2)%2)*LLR_lam_1[P.at(P.size()-1).at(5*Beta+3)] + pow(-1, u0)*LLR_lam_1[P.at(P.size()-1).at(5*Beta+4)];
					}
				}
				//printf("llr_%d[%d]=%lf\n", lam, Beta, LLR_lam[Beta]);
			}
		}
	}
}

void recursivelyUpdateC(int lam, int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<int> &kernel, vector<vector<int>> &P, int list_size){
	int psi;
	int N_p = 1;
	REP(i, kernel.size()-lam){
		N_p *= kernel.at(kernel.size() - 1 - i);
	}

	if(kernel.at(lam-1) == 2){
		psi = phi/2;

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			
			int *c_lam   = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int *c_lam_1 = getArrayPointer_C(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					c_lam_1[kernel.at(lam-2)*(2*Beta)   + (psi%kernel.at(lam-2))] = ( c_lam[2*Beta] + c_lam[2*Beta+1] )%2;
					c_lam_1[kernel.at(lam-2)*(2*Beta+1) + (psi%kernel.at(lam-2))] = c_lam[2 * Beta + 1];
				}
				else{
					c_lam_1[P.at(P.size()-1).at(2*Beta)]   = ( c_lam[2*Beta] + c_lam[2*Beta+1] )%2;
					c_lam_1[P.at(P.size()-1).at(2*Beta+1)] = c_lam[2 * Beta + 1];
				}
			}
		}
	}
	else if(kernel.at(lam-1) == 3){
		psi = phi/3;

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			
			int *c_lam   = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int *c_lam_1 = getArrayPointer_C(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					c_lam_1[kernel.at(lam-2)*(3*Beta)   + (psi%kernel.at(lam-2))] = ( c_lam[3*Beta] + c_lam[3*Beta+1] ) % 2;
					c_lam_1[kernel.at(lam-2)*(3*Beta+1) + (psi%kernel.at(lam-2))] = ( c_lam[3*Beta] + c_lam[3*Beta+2] ) % 2;
					c_lam_1[kernel.at(lam-2)*(3*Beta+2) + (psi%kernel.at(lam-2))] = ( c_lam[3*Beta] + c_lam[3*Beta+1] + c_lam[3*Beta+2] ) % 2;
				}
				else{
					c_lam_1[P.at(P.size()-1).at(3*Beta)]   = ( c_lam[3*Beta] + c_lam[3*Beta+1] ) % 2;
					c_lam_1[P.at(P.size()-1).at(3*Beta+1)] = ( c_lam[3*Beta] + c_lam[3*Beta+2] ) % 2;
					c_lam_1[P.at(P.size()-1).at(3*Beta+2)] = ( c_lam[3*Beta] + c_lam[3*Beta+1] + c_lam[3*Beta+2] ) % 2;
				}
			}
		}
	}
	else if(kernel.at(lam-1) == 5){
		psi = phi/5;

		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			
			int *c_lam   = getArrayPointer_C(lam, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			int *c_lam_1 = getArrayPointer_C(lam-1, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

			REP( Beta, N_p ){
				if(lam != 1){
					c_lam_1[kernel.at(lam-2)*(5*Beta)   + (psi%kernel.at(lam-2))] = ( c_lam[5*Beta] + c_lam[5*Beta+1] + c_lam[5*Beta+2] + c_lam[5*Beta+3] ) % 2;
					c_lam_1[kernel.at(lam-2)*(5*Beta+1) + (psi%kernel.at(lam-2))] = ( c_lam[5*Beta] + c_lam[5*Beta+3] ) % 2;
					c_lam_1[kernel.at(lam-2)*(5*Beta+2) + (psi%kernel.at(lam-2))] = ( c_lam[5*Beta] + c_lam[5*Beta+3] + c_lam[5*Beta+4] ) % 2;
					c_lam_1[kernel.at(lam-2)*(5*Beta+3) + (psi%kernel.at(lam-2))] = ( c_lam[5*Beta] + c_lam[5*Beta+2] + c_lam[5*Beta+4] ) % 2;
					c_lam_1[kernel.at(lam-2)*(5*Beta+4) + (psi%kernel.at(lam-2))] = ( c_lam[5*Beta] + c_lam[5*Beta+4] ) % 2;
				}
				else{
					c_lam_1[P.at(P.size()-1).at(5*Beta)]   = ( c_lam[5*Beta] + c_lam[5*Beta+1] + c_lam[5*Beta+2] + c_lam[5*Beta+3] ) % 2;
					c_lam_1[P.at(P.size()-1).at(5*Beta+1)] = ( c_lam[5*Beta] + c_lam[5*Beta+3] ) % 2;
					c_lam_1[P.at(P.size()-1).at(5*Beta+2)] = ( c_lam[5*Beta] + c_lam[5*Beta+3] + c_lam[5*Beta+4] ) % 2;
					c_lam_1[P.at(P.size()-1).at(5*Beta+3)] = ( c_lam[5*Beta] + c_lam[5*Beta+2] + c_lam[5*Beta+4] ) % 2;
					c_lam_1[P.at(P.size()-1).at(5*Beta+4)] = ( c_lam[5*Beta] + c_lam[5*Beta+4] ) % 2;
				}
			}
		}
	}

	if(psi%2==1 && kernel.at(lam-2)==2 || psi%3==2 && kernel.at(lam-2)==3 || psi%5==4 && kernel.at(lam-2)==5){
		recursivelyUpdateC( (int)(lam-1), psi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size);
	}
}

void continuePaths_UnfrozenBit(int phi, vector<bool> &activePath, vector<stack<int>> &inactiveArrayIndices, vector<vector<int>> &pathIndexToArrayIndex, vector<vector<int>> &arrayReferenceCount, vector<vector<double *>> &arrayPointer_LLR, vector<vector<int *>> &arrayPointer_C, vector<double> &pathMetric, stack<int> &inactivePathIndices, vector<int> &kernel, int list_size, vector<int *> &arrayPointer_info){
	vector<double> probForks(2*list_size);
	vector<double> prob;
	int i = 0;
	// populate probForks
	REP(l, list_size){
		if(activePath.at(l) == true){
			double *LLR_m = getArrayPointer_LLR(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			probForks.at(2*l)   = -( pathMetric.at(l) + log(1 + exp(-LLR_m[0])) );
			probForks.at(2*l+1) = -( pathMetric.at(l) + log(1 + exp( LLR_m[0])) );

			prob.push_back( probForks.at(2*l) );
			prob.push_back( probForks.at(2*l+1) );
			i++;
		}
		else{
			probForks.at(2*l)   = NAN;
			probForks.at(2*l+1) = NAN;
		}
	}
	int rho = min(2*i, list_size);
	vector<int> contForks(2*list_size);

	// START - implemention of .at(2) Algorithm13 line 14 -
	REP(l, 2*list_size){ contForks.at(l) = 0; }
	sort( prob.begin(), prob.end(), greater<double>() );

	double threshold = prob.at(rho-1);
	int num_paths_continued = 0;

	REP(l, 2*list_size){
		if(probForks.at(l) > threshold){
			contForks.at(l) = 1;
			num_paths_continued++;
		}
		if(num_paths_continued == rho){ break; }
	}

	if( num_paths_continued < rho ){
		REP(l, 2*list_size){
			if(probForks.at(l) == threshold){
				contForks.at(l) = 1;
				num_paths_continued++;
			}
			if(num_paths_continued == rho){ break; }
		}
	}
	// END - implemention of .at(2) Algorithm13 line 14 -

	// kill-off non-continuing paths
	REP(l, list_size){
		if(activePath.at(l) == 0){ continue; }
		if(contForks.at(2*l) == 0 && contForks.at(2*l+1) == 0){
			killPath(l, pathMetric, inactivePathIndices, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, kernel);
		}
	}

	// continue reveant paths, and duplicate if necessary
	REP(l, list_size){
		if( contForks.at(2*l)==0 && contForks.at(2*l+1)==0 ){ continue; }

		int * c_m = getArrayPointer_C(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

		if( contForks.at(2*l)==1 && contForks.at(2*l+1)==1 ){
			c_m[phi%kernel.at(kernel.size()-1)] = 0;
			int lp = clonePath(l, pathMetric, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, kernel);
			c_m = getArrayPointer_C(kernel.size(), lp, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			c_m[phi%kernel.at(kernel.size()-1)] = 1;
			copy(arrayPointer_info.at(l), arrayPointer_info.at(l) + phi,  arrayPointer_info.at(lp));
			arrayPointer_info.at(l)[phi] = 0;
			arrayPointer_info.at(lp)[phi] = 1;

			double *LLR_m = getArrayPointer_LLR(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			pathMetric.at(l)  += log( 1 + exp(-LLR_m[0]) );		// hat_u[i] = 0
			LLR_m = getArrayPointer_LLR(kernel.size(), lp, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
			pathMetric.at(lp) += log( 1 + exp( LLR_m[0]) );		// hat_u[i] = 1
		}
		else{
			if(contForks.at(2*l) == 1){
				c_m[phi%kernel.at(kernel.size()-1)] = 0;
				arrayPointer_info.at(l)[phi] = 0;
				double *LLR_m = getArrayPointer_LLR(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
				pathMetric.at(l) += log( 1 + exp(-LLR_m[0]) );
			}
			else{
				c_m[phi%kernel.at(kernel.size()-1)] = 1;
				arrayPointer_info.at(l)[phi] = 1;
				double *LLR_m = getArrayPointer_LLR(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
				pathMetric.at(l) += log( 1 + exp( LLR_m[0]) );
			}
		}
	}
}

void MK_SCL_Decoder(double variance, vector<double> &y, vector<int> &I, vector<int> &u_hat, vector<int> &kernel, vector<vector<int>> &P, int list_size, vector<double> &PMR){
	// initializeDataStructures()
	vector<double> pathMetric(list_size, 0.0);
	stack<int> inactivePathIndices;
	while(inactivePathIndices.size()){
		inactivePathIndices.pop();
	}
	vector<bool> activePath(list_size, 0);
	vector<int *> arrayPointer_info(list_size);
	vector<vector<double *>> arrayPointer_LLR(kernel.size()+1, vector<double *>(list_size));
	vector<vector<int *>> arrayPointer_C(kernel.size()+1, vector<int *>(list_size, 0));
	vector<vector<int>> pathIndexToArrayIndex(kernel.size()+1, vector<int>(list_size, 0));
	vector<stack<int>> inactiveArrayIndices;
	inactiveArrayIndices.resize(kernel.size()+1);
	REP(i, kernel.size()+1){
		while(inactiveArrayIndices[i].size()){ inactiveArrayIndices[i].pop(); }
	}
	vector<vector<int>> arrayReferenceCount(kernel.size()+1, vector<int>(list_size, 0));

	REP(s, list_size){
		arrayPointer_info.at(s) = new int[y.size()];
		REP(lam, kernel.size()+1){
			int N_p = 1;
			REP(i, kernel.size()-lam){
				N_p *= kernel.at(kernel.size() - 1 - i);
			}
			arrayPointer_LLR.at(lam).at(s) = new double[N_p]();
			arrayPointer_C.at(lam).at(s) = new int[(lam == 0 ? 1 : kernel.at(lam-1)) * N_p]();
			arrayReferenceCount.at(lam).at(s) = 0;
			inactiveArrayIndices.at(lam).push(s);
		}
	}
	REP(l, list_size){
		activePath.at(l) = false;
		inactivePathIndices.push(l);
	}

	int l = assignInitialPath(inactivePathIndices, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, kernel);
	double * LLR_0 = getArrayPointer_LLR(0, l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);

	//main loop of SCL Decoder
	REP(Beta, y.size()){    // LLR[0] initialization
		LLR_0[Beta] = 2 * y.at(Beta) / variance;
	}
	REP(phi, y.size()){
		recursivelyCalcLLR(kernel.size(), phi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size);
		if(any_of(I.begin(), I.end(), [&phi](int x) { return x == phi; })){		// information bit
			continuePaths_UnfrozenBit(phi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, pathMetric, inactivePathIndices, kernel, list_size, arrayPointer_info);
		}
		else{		// frozen bit
			REP(l, list_size){
				if(activePath.at(l) == false){ continue; }
				int * C_m = getArrayPointer_C(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
				double *LLR_m = getArrayPointer_LLR(kernel.size(), l, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel);
				pathMetric.at(l) += log( 1 + exp(-LLR_m[0]) );
				C_m[phi%kernel.at(kernel.size()-1)] = 0;
				arrayPointer_info.at(l)[phi] = 0;
			}
		}

		double min_PM = numeric_limits<double>::max();
		REP(l, list_size){
			if(activePath.at(l) == false){ continue; }
			if(pathMetric.at(l) < min_PM){ min_PM = pathMetric.at(l); }
		}
		PMR.at(phi) += *max_element(pathMetric.begin(), pathMetric.end()) - min_PM;		// add path metric range to PMR

		if(phi%kernel.at(kernel.size()-1) == kernel.at(kernel.size()-1)-1){
			recursivelyUpdateC(kernel.size(), phi, activePath, inactiveArrayIndices, pathIndexToArrayIndex, arrayReferenceCount, arrayPointer_LLR, arrayPointer_C, kernel, P, list_size);
		}
	}

	// Return the best codeword in the list
	int lp = 0;
	double p_llr = numeric_limits<double>::max();
	REP(l, list_size){
		if(activePath.at(l) == false){ continue; }
		
		if(pathMetric.at(l) < p_llr){
			p_llr = pathMetric.at(l);
			lp = l;
		}
	}

	int * C_0 = arrayPointer_info.at(lp);
	REP(Beta, u_hat.size()){
		u_hat.at(Beta) = C_0[Beta];
	}

	// delete data
	REP(s, list_size){
		delete arrayPointer_info.at(s);
		REP(lam, kernel.size()+1){
			delete arrayPointer_LLR.at(lam).at(s);
			delete arrayPointer_C.at(lam).at(s);
		}
	}
}