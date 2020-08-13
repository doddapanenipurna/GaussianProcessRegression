#ifndef GPR_H
#define GPR_H
#include "math.h"
#include "iostream"
#include "omp.h"

using namespace std;
class GPR{
	public:
		~GPR(){};
		void run_GPR(double* ftest, double** XY, double* f, int* itest, int* itrain, double l1, double l2, int itest_size, int itrain_size, double t, int n);	
		double kernel_comp(double x, double y, double u, double v, double l1, double l2);
		double** kernel(double** XY, double L1, double L2, int n);
		double** data_extract(double** K0, int* itrain, int n);
		void print_matrix(double** myMatrix, int size);
		double** transpose(double** myMatrix, int n);
		double* backward_substitution(double** A, double* b, int n);
		double* forward_substitution(double** A, double* b, int n);
		double** chol_fact(double** K, int n, double t);
		double** data_extract2(double** K0, int* itrain, int* itest, int ntrain, int ntest);
};
#endif
