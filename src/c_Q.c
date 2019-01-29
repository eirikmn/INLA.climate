#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
void c_Q(double *ii,double *jj,double *xx,int n,int m,double *ww,double *pp,double tau,double sx)
{
	
	for(int i=0;i<n*(m+1);i+=1){
		ii[i] = i/(m+1)+1;
		jj[i] = 1+(i%(m+1))*n +i/(m+1);
		
		if(i%(m+1)==0){
			xx[i] = tau;
		}else{
			xx[i] = -sx*sqrt(ww[i%(m+1)-1])*tau;
		}
	}
	int teller = n*(m+1);
	int k=0;
	
	for(int i=0;i<m*n;i+=1){
		ii[teller] = i+n+1;
		jj[teller] = i+n+1;
		k=i/n;
		xx[teller]= 1/(1-pp[k]*pp[k]) + sx*sx*ww[k]*tau;
		if(i%n!=0 && i%n!=n-1){
			xx[teller]=xx[teller]+pp[k]*pp[k]/(1-pp[k]*pp[k]);
		}
		teller += 1;
	}
	for(int k=0;k<m;k+=1){
		for(int i=0;i<(n-1);i+=1){
			ii[teller] = (k+1)*n+i+1;
			jj[teller] = (k+1)*n+i+2;
			xx[teller] = -pp[k]/(1-pp[k]*pp[k]);
			teller+=1;
		}
	}
	for(int bi=0;bi<m;bi+=1){
		for(int bj=bi+1;bj<m;bj+=1){
			for(int i=0;i<n;i+=1){
				ii[teller]= (bi+1)*n+i+1;
				jj[teller]= (bj+1)*n+i+1;
				xx[teller]= sx*sx*sqrt(ww[bi]*ww[bj])*tau;
				teller+=1;
			}
		}
	}
	
	
}
