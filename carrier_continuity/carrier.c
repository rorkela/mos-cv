#include "carrier.h"
#include "../main.h"
#define B(z) (z) / (exp((z)) - 1)
void carrier_continuity(double *V, double *Vprev, double *nprev, double *pprev, double *n, double *p, int N) {
  // Newton Rhapson used
  // Jac is Jacobian
  // Res  is the value of function from initial guess
  double *jac = calloc(3 * N, sizeof(double));
  double *res = malloc(N * sizeof(double));
  double *update = malloc(N*sizeof(double));
  int maxiter=20;
  int iter=0;
  do{
    computeJacobi_n(jac,mos.mu_n,V,p,N);
    residual_n(res, n, p, nprev, pprev, V, Vprev, mos.mu_n, N);
    for(int i=0;i<N;i++) res[i]=-res[i];
    thomas(jac,res,N,update);
    for(int i=0;i<N;i++) n[i]+=update[i];
  }while(iter++<maxiter);
  iter=0;
  do{
    computeJacobi_p(jac,mos.mu_p,V,n,N);
    residual_p(res, n, p, nprev, pprev, V, Vprev, mos.mu_p, N);
    for(int i=0;i<N;i++) res[i]=-res[i];
    thomas(jac,res,N,update);
    for(int i=0;i<N;i++) p[i]+=update[i];
  }while(iter++<maxiter);
  free(jac);
  free(res);
  free(update);
}
// J[i] is current density at i+0.5.
// WARNING: N+0.5 and -0.5 are not considered yet. I am assuming zero for now.
void compute_J(double *J, double *V, double *n, double u, int N) {
  for (int i = 0; i < N - 2; i++) {
    if(IN_OX(i))
      J[i]=0;
    else
      J[i] = (kB * mos.T / mos.dx) * u * B(V[i + 1] - V[i]) * (n[i + 1] - n[i] * exp(V[i + 1] - V[i]));
  }
  J[N-1]= (kB * mos.T / mos.dx) * u * 1 * ( n_teq - n[N-1] * 1);
  return;
}
// unchanged above part have to fix it if needed


// assuming the flow is V->N->P

// function that takes the matrices needed and throws back the jacobi of it
// rn this one is only for N but can be simply copied for P 
// INCOMPLETE and dont run
// u is the mobility
void computeJacobi_n(double *Jac,double u,double *V,double *p, int N){
  // so fill the jacobi with the required derivatives
  // jacobi is Nx3 as tridiagonal and the terms will come using the other matrices

  // WARNING: the term del Ri/del ni doesnt have the derivative of the generation and recombination rate wrt Ni

  double c=(kB*mos.T/(2*mos.dx*mos.dx));
  // just put cu/q as mogger_constant
  // not valid constant for a metal oxide
  // different mogger constant for p
  double mogger_constant= -(c*u/q);

  double delta_T=sim.dt;

  for(int i=0;i<N-1;i++){

    if(IN_OX(i)==1){
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;
    }

    else{
    Jac[3*i]= -(mogger_constant*(B(V[i]-V[i-1]) * exp(V[i]-V[i-1])));
    // DRi/Dni
    // WARNING: recombination derivative term left
    Jac[3*i+1]= (1/delta_T) + (mogger_constant*( (B(V[i+1]-V[i])*exp(V[i+1]-V[i]) ) + B(V[i]-V[i-1]) )) + mos.C_Rr*p[i]/2 ;
    // DRi/Dni+1
    Jac[3*i+2]= -(mogger_constant*(B(V[i+1]-V[i])));
    }

  }
  // boundary conditions 

  Jac[0]=0;//the first derivative is always zero as doesnt exist for a tri diagonal matrix

  //assumption that it is n doped  
  //WARNING: Ohmic coundition hardcoded
  Jac[3*(N-1)+1]=(1/delta_T) + (mogger_constant*(1+B(V[(N-1)]-V[(N-1)-1])))+ mos.C_Rr*p[N-1]/2 ;//ohmic boundary condition
  Jac[3*(N-1)+2]=0;
  Jac[3*(N-1)]= -(mogger_constant*(B(V[(N-1)]-V[(N-1)-1]) * exp(V[(N-1)]-V[(N-1)-1])));
}

void computeJacobi_p(double *Jac,double u,double *V,double *n, int N){
  // so fill the jacobi with the required derivatives
  // jacobi is Nx3 as tridiagonal and the terms will come using the other matrices

  // WARNING: the term del Ri/del ni doesnt have the derivative of the generation and recombination rate wrt Ni

  double c=(kB*mos.T/(2*mos.dx*mos.dx));
  // just put cu/q as mogger_constant
  // not valid constant for a metal oxide
  // different mogger constant for p
  double mogger_constant= (c*u/q);

  double delta_T=sim.dt;

  for(int i=0;i<N-1;i++){

    if(IN_OX(i)==1){
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;
    }

    else{
    Jac[3*i]= -(mogger_constant*(B(V[i]-V[i-1]) * exp(V[i]-V[i-1])));
    // DRi/Dni
    // WARNING: recombination derivative term left
    Jac[3*i+1]= (1/delta_T) + (mogger_constant*( (B(V[i+1]-V[i])*exp(V[i+1]-V[i]) ) + B(V[i]-V[i-1]) )) + mos.C_Rr*n[i]/2 ;
    // DRi/Dni+1
    Jac[3*i+2]= -(mogger_constant*(B(V[i+1]-V[i])));
    }

  }
  // boundary conditions 

  Jac[0]=0;//the first derivative is always zero as doesnt exist for a tri diagonal matrix

  //assumption that it is n doped  
  //WARNING: Ohmic coundition hardcoded
  Jac[3*(N-1)+1]=(1/delta_T) + (mogger_constant*(1+B(V[(N-1)]-V[(N-1)-1])))+ mos.C_Rr*n[N-1]/2 ;//ohmic boundary condition
  Jac[3*(N-1)+2]=0;
  Jac[3*(N-1)]= -(mogger_constant*(B(V[(N-1)]-V[(N-1)-1]) * exp(V[(N-1)]-V[(N-1)-1])));
}

void residual_n(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N){
  // function to return the values of residual
  // used crank nicolson scheme for this (can be changed if not working)
  // we need current density for this
  double *J = malloc(N * sizeof(double));
  compute_J(J,V,n,u,N);
  double *Jprev = malloc(N * sizeof(double));
  compute_J(Jprev,Vprev,nprev,u,N);

  for(int i=0;i<N;i++){
    // inside oxide and metal assumed no change in n and p and took values as zero
    if(IN_OX(i)){
      res[i]=0;
    }
    
    else{
      // sign of electron is handled
      // WARNING: the scheme changes for the function definition of compute_J
      res[i]=(n[i]-nprev[i])/(sim.dt) + ((J[i]-J[i-1]+Jprev[i]-Jprev[i-1])/(2*q*mos.dx)) -(mos.Gr - mos.C_Rr*(n[i]*p[i]/2+nprev[i]*pprev[i]/2-mos.ni*mos.ni));
    }
  }
  free(J);
  free(Jprev);
}

void residual_p(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N){
  // function to return the values of residual
  // used crank nicolson scheme for this (can be changed if not working)
  // we need current density for this
  double *J = malloc(N * sizeof(double));
  compute_J(J,V,p,u,N);
  double *Jprev = malloc(N * sizeof(double));
  compute_J(Jprev,Vprev,pprev,u,N);

  for(int i=0;i<N;i++){
    // inside oxide and metal assumed no change in n and p and took values as zero
    if(IN_OX(i)){
      res[i]=0;
    }
    
    else{
      // sign of electron is handled
      // WARNING: the scheme changes for the function definition of compute_J
      res[i]=(n[i]-nprev[i])/(sim.dt) + ((J[i]-J[i-1]+Jprev[i]-Jprev[i-1])/(2*q*mos.dx)) -(mos.Gr - mos.C_Rr*(n[i]*p[i]/2+nprev[i]*pprev[i]/2-mos.ni*mos.ni));
    }
  }
  free(J);
  free(Jprev);
}


