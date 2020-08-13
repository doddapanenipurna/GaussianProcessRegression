#ifndef GENERATORS_H
#define GENERATORS_H

#include <iostream>
#include "math.h"
using namespace std;

class Generators{
	public:
		~Generators(){};
		double** generate_XY(int m);
		void print_XY(double** myMatrix, int size);
		double* generate_Kernel_1d(double** XY, double* Y, double L1, double L2, int m);
		double* generate_F(double** XY, int m);
		void print_F(double* myMatrix, int size);
		double* modify_f(double** XY, double* f, double* kernel, double xval, double yval, int m);
		int* random_perm(int n);
		int* take_slice(int* arr, int start, int end);
		double* generate_L(int m, double min, double max, double step);
		double** generate_MSE(int n);
};
#endif
