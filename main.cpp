#include <iostream>
#include "generators.h"
#include "gpr.h"
using namespace std;

int main(int argc, char* argv[]){
	int m = atoi(argv[1]);
	int n = m*m;
	int ntest = round(0.1*n);
	int ntrain = n - ntest;
	//decalre new generator
	//instantiate neccesary data vectors and matrices
	Generators* myGenerators = new Generators();
	GPR* myGPR = new GPR();
	double** XY = myGenerators->generate_XY(m);
	double* F = myGenerators -> generate_F(XY, m);
	double temp[2] = {0.25,0.25};
	double* Y = temp;
	double* kernel = myGenerators -> generate_Kernel_1d(XY,Y, 2.0/m,2.0/m,m); 
	F = myGenerators -> modify_f(XY, F, kernel, 0.2, 0.1, n);
	int* rperm = myGenerators -> random_perm(n);
	int* itest = myGenerators -> take_slice(rperm, 0, ntest);
	int* itrain = myGenerators -> take_slice(rperm, ntest+1, n);
	
	double Tparam = 0.5;
	double* Lparam = myGenerators -> generate_L(m, 0.25, 10.0, 0.5);
	int L_size = (10.0-0.25)/0.5;
	L_size--;
	//L_size++;
	/*
	for(int x = 0; x<L_size; x++){
		cout<<Lparam[x]<<"\n";
	}
	for(int x = 0; x< ntrain; x++){
		cout<<itrain[x]<<"\n";
	}*/
	double start = 0.0;
        double end = 0.0;	
	start = omp_get_wtime();
	
        double* ftest = new double[ntrain];
        double* error = new double[ntrain];
        double** MSE = myGenerators -> generate_MSE(L_size);
	for(int x = 0; x< L_size; x++){
		for(int y = 0; y< L_size; y++){
				myGPR-> run_GPR(ftest, XY, F, itest, itrain, Lparam[x], Lparam[y], ntest, ntrain,Tparam, n);
			double MSE_val = 0.0;
			for(int z = 0; z<ntest; z++){
				error[z] = F[itest[z]]-ftest[z];
				error[z] = error[z]*error[z];
				MSE_val += error[z];
			}
			MSE[x][y] = MSE_val/ntest;
			cout<<"l1: "<<Lparam[x]<<", l2: "<<Lparam[y]<<", MSE: ";
			cout<<MSE[x][y]<<"\n";
		}
	}
	end = omp_get_wtime();
	cout<<"Final Time: "<<end-start<<"\n";
}
