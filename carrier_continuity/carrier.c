#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../parameter_fetch/parameter_fetch.h"
#include "carrier.h"
#define B(z) (z)/(exp((z))-1)
void carrier_continuity(double * V, double * nprev, double * pprev, double * n, double * p, int N)
{
    // Newton Rhapson used
    // Jac is Jacobian
    // F  is the value of function from initial guess 
    double *Jac=calloc(N*N,sizeof(double));
    double *F=malloc(N*sizeof(double));
    // initializing current
    double *Jpprev=malloc((N-1)*sizeof(double));
    double *Jp=malloc((N-1)*sizeof(double));
    double *Jnprev=malloc((N-1)*sizeof(double));
    double *Jn=malloc((N-1)*sizeof(double));
    compute_J(Jpprev,V,pprev,mos.mu_p,N);
    compute_J(Jp,V,p,mos.mu_p,N);
    compute_J(Jnprev,V,nprev,mos.mu_n,N);
    compute_J(Jn,V,n,mos.mu_n,N);
    free(Jpprev);
    free(Jp);
    free(Jnprev);
    free(Jn);
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
