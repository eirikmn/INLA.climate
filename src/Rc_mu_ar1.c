
void c_mu_ar1(double*,double *,int,int,double *, double *,double,double);
void Rc_mu_ar1(double *x,double *forcing, int *n,int *m, double *w, double *lambda,double *sf, double *F0)
{
		
   c_mu_ar1(x,forcing,*n,*m,w,lambda,*sf,*F0);
   
   
}
