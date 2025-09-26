#include "../main.h"
#define B(z) (z) / (exp((z)) - 1)
void carrier_continuity(double *V, double *Vprev, double *nprev, double *pprev, double *n, double *p, int N) {
  // Newton Rhapson used
  // Jac is Jacobian
  // F  is the value of function from initial guess
  double *Jac_n = calloc(3 * N, sizeof(double));
  double *Jac_p = calloc(3 * N, sizeof(double));
  double *Res = malloc(N * sizeof(double));
  // Solving for n
  //compute_F(F, V, Vprev, nprev, pprev, n, p, N, 1);
  //compute_F(F, V, Vprev, nprev, pprev, n, p, N, 0);
  free(Jac_n);
  free(Jac_p);
  free(Res);
}
// J[i] is current density at i+0.5.
// WARNING: N+0.5 and -0.5 are not considered yet. I am assuming zero for now.
void compute_J(double *J, double *V, double *n, double u, int N) {
  for (int i = 0; i < N - 1; i++) {
    J[i] = (kB * mos.T / mos.dx) * u * B(V[i + 1] - V[i]) * (n[i + 1] - n[i] * exp(V[i + 1] - V[i]));
  }
  return;
}
// unchanged above part have to fix it if needed
// assuming the flow is V->N->P

// function that takes the matrices needed and throws back the jacobi of it
// rn this one is only for N but can be simply copied for P 
// INCOMPLETE and dont run
void computeJacobi_n(double *Jac,double* n,double* nprev,double* p,double* pprev,double u,double *V, double *Vprev, int N){
  // so fill the jacobi with the required derivatives
  // jacobi is Nx3 as tridiagonal and the terms will come using the other matrices

  // NOTE: three things havent been considered rn
  // 1) the boundary conditions
  // 2) the respective jacobi entries for the metal oxide region has to be implemented
  // 3) the term del Ri/del ni doesnt have the derivative of the generation and recombination rate wrt Ni

  double c=(kB*mos.T/mos.dx);
  // just put cu/q as mogger_constant
  double mogger_constant= (c*u/q);

  double delta_T=;

  for(int i=0;i<N;i++){

    // DRi/Dni-1
    Jac[3*i]= -(mogger_constant*(B(V[i]-V[i-1]) * exp(V[i]-V[i-1])));
    // DRi/Dni
    Jac[3*i+1]= (1/delta_T) + (mogger_constant*( (B(V[i+1]-V[i])*exp(V[i+1]-V[i]) ) + B(V[i]-V[i-1]) )) ;
    // DRi/Dni+1
    Jac[3*i+2]= -(mogger_constant*(B(V[i+1]-V[i])));

    // boundary conditions and then the metal oxide considerations
    Jac[0]=0;//the first derivative is always zero as doesnt exist for a tri diagonal matrix
    Jac[3*(N-1)+2]=;//ohmic boundary condition


  }
}


