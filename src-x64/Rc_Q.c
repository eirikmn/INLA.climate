void c_Q(double*,double*,double*,int,int,double*,double*,double,double);
void Rc_Q(double* ii,double* jj,double* xx,int* n,int* m,double* ww,double* pp,double* tau, double* sx){
    c_Q(ii,jj,xx,*n,*m,ww,pp,*tau,*sx);
}