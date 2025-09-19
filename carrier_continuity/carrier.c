#include "../main.h"
#define B(z) (z)/(exp((z))-1)
void carrier_continuity(double * V, double * Vprev, double * nprev, double * pprev, double * n, double * p, int N)
{
    // Newton Rhapson used
    // Jac is Jacobian
    // F  is the value of function from initial guess 
    double *Jac=calloc(3*N,sizeof(double));
    double *F=malloc(N*sizeof(double));
    //Solving for n
    compute_F(F,V,Vprev,nprev,pprev,n,p,N,1);
    compute_F(F,V,Vprev,nprev,pprev,n,p,N,0);
    free(Jac);
    free(F);
}
// J[i] is current density at i+0.5.
// WARNING: N+0.5 and -0.5 are not considered yet. I am assuming zero for now.
void compute_J(double *J, double *V, double *n,double u, int N)
{
    for(int i=0;i<N-1;i++)
    {
        J[i]=(kB*mos.T/mos.dx)*u*B(V[i+1]-V[i])*(n[i+1]-n[i]*exp(V[i+1]-V[i]));
    }
    return;
}
// INCOMPLETE
void compute_F(double *F, double * V, double * Vprev, double *nprev, double *pprev, double * n, double * p, int N, int is_n)
{
    double u=mos.mu_p;
    double qfact=1;
    double F1=0;
    double F2=0;
    double J0=0;
    double J1=0;
    if(is_n==1)
    {
        u=mos.mu_n;
        qfact=-1;
    }
    for(int i=0;i<N;i++)
    {
        if(IN_OX(i))
        {
            F[i]=0;
        }
        else
        {
            J0=(kB*mos.T/mos.dx)*u*B(V[i+1]-V[i])*(n[i+1]-n[i]*exp(V[i+1]-V[i]));
        }
    }
}
