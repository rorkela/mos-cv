#include "../main.h"
#define MAX_ITER 20
double solve_c(struct signal Vin)
{
    int N=mos.nz;
    double drichlet_factor=kB*mos.T*log(mos.Nc/mos.Nd); // WARNING: Hardcoded for n doped
    double delta=0;
    int iter=0;
    // Initializing arrays for n p V and for previous time instant
    double *n=malloc(N*sizeof(double));
    double *p=malloc(N*sizeof(double));
    double *V=malloc(N*sizeof(double));
    double *n_prev_t=malloc(N*sizeof(double));
    double *p_prev_t=malloc(N*sizeof(double));
    double *V_prev_t=malloc(N*sizeof(double));

    //Defining charge density
    double *charge_density=malloc(N*sizeof(double));

    // Starting from thermal equilibrium conditions at t=0;
    for(int i=0;i<N;i++)
    {
        if(IN_OX(i))
        {
            n[i]=0;
            p[i]=0;
        }
        // WARNING: Assuming this is n doped. Change here if needed. Hardcoded for now
        else{
        n[i]=mos.Nd; 
        p[i]=mos.ni*mos.ni/mos.Nd;
        }
        V[i]=Vin.bias*(1-(double)i/(N-1)) + drichlet_factor;
    }
    // NOTE: DEBUG
    plotxy(sim.x, V, N);
    //Storing copy for t=0.
    copy_arr(n,n_prev_t,N);
    copy_arr(p,p_prev_t,N);
    copy_arr(V,V_prev_t,N);

    do{
        poisson(V,n,p,V[0],V[N-1]);
        carrier_continuity(V,V_prev_t,n_prev_t,p_prev_t,n,p,N);
        
        // Logic for computing delta
        delta=0;
        for(int i=0;i<N;i++)
        {
            compute_delta(&delta,V[i],V_prev_t[i]);
            compute_delta(&delta,n[i],n_prev_t[i]);
            compute_delta(&delta,p[i],p_prev_t[i]);
        }
        if(iter++>=MAX_ITER) break;
    }while(delta>=1e-6);
        plotxy(sim.x, V, N);

    solve_charge_density(charge_density,n,p,N);
    printf("Solved");
    free(n);
    free(p);
    free(V);
    return 0;
}

void solve_charge_density(double *charge_density, double *n, double *p, int N) //To solve for charge density
{
    for(int i=0;i<N;i++)
    {
        charge_density[i]=q*(p[i]-n[i]+mos.Nd-mos.Na);
    }
}

void copy_arr(double * source , double * target , int N)
{
    for(int i=0;i<N;i++) target[i]=source[i];
}
void compute_delta(double *delta, double val, double valprev)
{
    
}
