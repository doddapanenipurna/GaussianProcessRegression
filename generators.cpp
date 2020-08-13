#include "generators.h"

double** Generators::generate_XY(int m){
        int n = m*m;
        double h = 1.0/(m+1.0);
        double** XY = new double*[n];
        for(int x =0; x< n; x++){
                        XY[x]=new double[2];
        }
        int idx = 0;
        for(int i=0; i<m;i++){
                for(int j=0; j<m;j++){
                        idx++;
                        XY[idx-1][0] = (i+1)*h;
                        XY[idx-1][1] = (j+1)*h;
                }
        }
        return(XY);
}


void Generators::print_XY(double** myMatrix, int size){
        for(int x = 0; x<size*size; x++){
                cout<<myMatrix[x][0]<<",";
                cout<<myMatrix[x][1]<<"\n";
        }
}

double* Generators::generate_Kernel_1d(double** XY, double* Y, double L1, double L2, int m){
	double pi_const = 1/sqrt(2*M_PI);
	double* kernel = new double[m*m];
	for(int x =0; x< m*m; x++){	
		double xs = 0.0;
		double ys = 0.0;
		double my_exp = 0.0;
		xs = XY[x][0]-Y[0];
		xs = (xs*xs)/(2*(L1*L1));
		ys = XY[x][1]-Y[1];
		ys = (ys*ys)/(2*(L2*L2));
		my_exp = -1*(xs+ys);
		kernel[x] = pi_const * exp(my_exp);
	}
	return(kernel);	
}


double* Generators::generate_F(double** XY, int m){
        int n = m*m;
        double* f = new double[n];
        for(int x=0; x< n; x++){
                f[x] = 0.02*(((double)rand()/RAND_MAX)-0.5);
        }
        return(f);
}

void Generators::print_F(double* myMatrix, int size){
        cout<<"f = [";
        for(int x = 0; x<size; x++){
                if(x == size -1){
                        cout<<myMatrix[x]<<"\n";
                }else{
                 cout<<myMatrix[x]<<";\n";
                }
        }
        cout<<"]";
}

double*  Generators::modify_f(double** XY, double* f, double* kernel, double xval, double yval, int n){
	for(int x = 0; x< n; x++){
		f[x] += kernel[x];
		f[x] += XY[x][0]*xval + XY[x][1]*yval;
		
	}
	return(f);
}

int* Generators::random_perm(int n){
	int* myRandArr = new int[n];
	for(int x =0; x<n; x++){
		myRandArr[x] = x;
	}
	for(int i =0 ; i<n; i++){
		int j,t;
		j = rand() % (n-i) + i;
		t = myRandArr[j];
		myRandArr[j] = myRandArr[i];
		myRandArr[i] = t;	
	}
	return(myRandArr);
}

int* Generators::take_slice(int* arr, int start, int end){
	int* sliceArr = new int[end-start];
	int y =0;
	for(int x = start; x< end; x++){
		sliceArr[y] = arr[x];
		y++;
	}
	return(sliceArr);
}
	
double* Generators::generate_L(int m, double min, double max, double step){
	double size1 = round((max+1 - min)/step);
	int size = (int)size1;
	double* L_param = new double[size+1];
	int y = 0;
	for(double x = min; x<=max+1; x+= step){
		L_param[y] = x/m;
		y++;	
	}
	return(L_param);
}

double** Generators::generate_MSE(int n){
	double** MSE = new double*[n]; 
        for(int x =0; x< n; x++){
                        MSE[x]=new double[n];
        }
	for(int x = 0; x<n; x++){
		for(int y=0; y<n; y++){
			MSE[x][y] = 0;
		}
	}
	return(MSE);
}
