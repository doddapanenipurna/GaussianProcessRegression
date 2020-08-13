#include "gpr.h"
#include "omp.h"

void GPR::print_matrix(double** myMatrix, int size){
        for(int x =0; x<size; x++){
                for(int y =0; y<size;y++){
                        cout<<myMatrix[x][y]<<" ";
                }
        cout<<"\n";
        }
        cout<<"___________________________________\n";
}

double** GPR::kernel(double** XY, double l1, double l2, int n){
	double x,u,y,v =0;
	double** K0 = new double*[n];
	for(int x =0; x< n; x++){
		K0[x]=new double[n];
        }
	for(int i = 0; i<n; i++){
		x= XY[i][0];
		y= XY[i][1];
		for(int j = 0; j<n; j++){
			u = XY[j][0];
			v = XY[j][1];
			K0[i][j] = this->kernel_comp(x,y,u,v,l1,l2);
		}
	}
	return(K0);
}

double GPR::kernel_comp(double x, double y, double u, double v, double l1, double l2){
	double pi_const = 1.0/sqrt(2.0*M_PI);
	double my_exp = ((x-u)*(x-u))/(2.0*(l1*l1));
	my_exp += ((y-v)*(y-v))/(2.0*(l2*l2));
	my_exp = -1.0*my_exp;
	double result = pi_const * exp(my_exp);
	return(result);
}

double** GPR::data_extract(double** K0, int* itrain, int n){	
	double** K = new double*[n];
	for(int x = 0; x<n; x++){
		K[x] = new double[n];
	}
	int xarr = 0;
	int yarr =0;
	for(int x = 0; x<n; x++){
		xarr = itrain[x];
		for(int y = 0; y<n; y++){
			yarr = itrain[y];
			K[x][y] = K0[xarr][yarr];
		}
	}
	return(K);
}


double** GPR::data_extract2(double** K0, int* itrain, int* itest, int ntrain, int ntest){
	double** K = new double*[ntrain];
	for(int x = 0; x<ntrain; x++){
		K[x] = new double[ntest];
	}
	int xarr = 0;
	int yarr = 0;
	
	for(int x = 0; x<ntrain; x++){
		xarr = itrain[x];
		for(int y = 0; y<ntest; y++){
			yarr = itest[y];
			K[x][y] = K0[xarr][yarr];
		}
	}
	return(K);
	
}
void GPR::run_GPR(double* ftest, double** XY, double* f, int* itest, int* itrain, double l1, double l2, int itest_size, int itrain_size, double t, int n){

	double** K0 = this -> kernel(XY,l1,l2,n);
	double** train = this -> data_extract(K0, itrain, itrain_size);
	double* f_train = new double[itrain_size];
	for(int x = 0; x< itrain_size; x++){
		f_train[x] = f[itrain[x]];
	}
	double** L = chol_fact(train,itrain_size,t);
	double* F = forward_substitution(L,f_train,itrain_size);
	L = transpose(L,itrain_size);
	F = backward_substitution(L,F,itrain_size);
	double** k = this -> data_extract2(K0, itrain, itest, itrain_size, itest_size); 

	for(int x = 0; x<itrain_size;x++){
		ftest[x] = 0;
		for(int y = 0; y<itest_size;y++){
                	ftest[x] += k[x][y]*F[x];
		}
        }
}

double** GPR::chol_fact(double** K, int n, double t){
        for(int x =0; x<n; x++){
                K[x][x] = K[x][x]+t;
        }

        int a, b ,c;
	omp_lock_t lock;
	omp_init_lock(&lock);

        for(b= 0; b<n;b++){
                for(a=0;a<b;a++){
                        K[a][b] = 0;
                }
		#pragma omp parallel for shared(K),private(c)
                for(c = 0; c<a;c++){
                        K[b][b] = K[b][b] - K[b][c]*K[b][c];
                }
		#pragma omp single
                K[a][a] = sqrt(K[b][b]);
                #pragma omp parallel for shared(K), private(a,c)
		for(a = b+1; a<n;a++){
                        for(c=0;c<b;c++){
                                K[a][b] = K[a][b] - K[a][c]*K[b][c];
                        }
                        K[a][b] = K[a][b]/K[b][b];
                }
		omp_destroy_lock(&lock);

        }
        return(K);
}
double* GPR::forward_substitution(double** A, double* b, int n){
        double* x = new double[n];
        for(int y = 0; y<n;y++){
                x[y] = 0.0;
        }
        int i,j;
        double sum;
	#pragma omp parallel for private(i,j,sum)
        for(i =0; i<n;i++){
                sum = 0.0;
                for(j = 0; j<i;j++){
                        sum += A[i][j]*x[j];
                }
                x[i]=(b[i]-sum)/A[i][i];
        }
        return(x);
}


double* GPR::backward_substitution(double** A, double* b, int n){
        double* x = new double[n];
        x[n-1] = b[n-1]/A[n-1][n-1];
        int  i, j;
        double sum;
        #pragma omp parallel for private(i,j,sum) 
	for(i = n-2; i>-1; i--){
                sum = 0.0;
                for(j = i+1; j<n; j++){
                        sum = sum + A[i][j]*x[j];
                }
                sum = b[i]- sum;
                x[i] = sum/A[i][i];
        }
        return(x);
}

double** GPR::transpose(double** myMatrix, int n){
        double** newMatrix = new double*[n];
        for(int x =0; x< n; x++){
                        newMatrix[x]=new double[n];
        }
        for(int x = 0; x<n; x++){
                for(int y= 0; y<n; y++){
                        newMatrix[x][y] = myMatrix[y][x];
                }
        }
        return(newMatrix);
}

