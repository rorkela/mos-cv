#include "solve_c.h"
double solve_c(struct signal Vin)
{
    // Generating initial value for n and p
    int N=mos.nz;
    double *n=malloc(N*sizeof(double));
    double *p=malloc(N*sizeof(double));
    double temp=mos.ni*mos.ni/mos.Na;
    for(int i=0;i<N;i++)
    {
        p[i]=mos.Na;
        n[i]=temp;
    }
    //Defining charge density
    double *charge_density=malloc(N*sizeof(double));
    
    solve_charge_density(charge_density,n,p,N);
    printf("Solved");
    return 0;
}

void solve_charge_density(double *charge_density, double *n, double *p, int N) //To solve for charge density
{
    for(int i=0;i<N;i++)
    {
        charge_density[i]=q*(p[i]-n[i]+mos.Nd-mos.Na);
    }
}
