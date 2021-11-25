#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

void c_mu_ar1(double *means,double *forcing,int n,int m, double *w, double *lambda,double sf,double F0)
{
	
	
	
	double* z = (double*)R_alloc(n,sizeof(double));
    double* structure = (double*)R_alloc(n,sizeof(double));

	
	for(int i=0;i<n;i=i+1){
        z[i] = sf*(forcing[i]+F0);
		structure[i] = 0;
		
		for(int k=0;k<m;k=k+1){
			structure[i] += w[k]*pow(M_E,lambda[k]*(0.5+(double)i));
		}
		
	}
    
    for(int i=0;i<n;i=i+1){
        for(int j=0;j<=i;j=j+1){
			means[i] +=    z[j]*structure[i-j];
		}
    }
	
	//R_free(z);
	//R_free(structure);
	//free(z);
	//free(structure);
	
	
}
