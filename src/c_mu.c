#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

void c_mu(double *means,double *forcing,int n,double H,double sf,double shift)
{
	
	
	
	double* z = (double*)R_alloc(n,sizeof(double));
    double* structure = (double*)R_alloc(n,sizeof(double));

	
	for(int i=0;i<n;i=i+1){
        z[i] = sf*(forcing[i]+shift);
    }
	
    
    for(int i=0;i<n;i=i+1){
        structure[i]=pow(0.5+(double)i,H-1.5);
        for(int j=0;j<=i;j=j+1){

			means[i] +=    z[j]*structure[i-j];
        }
		
    }
	
	//R_free(z);
	//R_free(structure);
	//free(z);
	//free(structure);
	
	
}
