#include "carrier.h"
#include "../main.h"
#define B(z) (z==0.0?1:(z) / expm1((z)))
void carrier_continuity(double *V, double *Vprev, double *nprev, double *pprev, double *n, double *p, int N) {
  // Newton Rhapson used
  // Jac is Jacobian
  // Res  is the value of function from initial guess
  double *jac = malloc(3*N*sizeof(double));
  double *res = malloc(N*sizeof(double));
  double *update = malloc(N*sizeof(double));
  double *Vnorm =  malloc(N*sizeof(double));
  double *Vprevnorm = malloc(N*sizeof(double));
  for(int i=0;i<N;i++) Vnorm[i]=V[i]/(kB*mos.T/q);
  for(int i=0;i<N;i++) Vprevnorm[i]=Vprev[i]/(kB*mos.T/q);
  int maxiter=10;
  int iter=0;
  computeJacobi_n(jac,mos.mu_n,Vnorm,p,N);
  do{
    residual_n(res, n, p, nprev, pprev, Vnorm, Vprevnorm, mos.mu_n, N);
    for(int i=0;i<N;i++) res[i]=-res[i];
    thomas(jac,res,N,update);
    for(int i=0;i<N;i++) n[i]+=0.2*update[i];
  }while (iter++<maxiter);
  iter=0;
  computeJacobi_p(jac,mos.mu_p,Vnorm,n,N);
  do{
    residual_p(res, n, p, nprev, pprev, Vnorm, Vprevnorm, mos.mu_p, N);
    for(int i=0;i<N;i++) res[i]=-res[i];
    thomas(jac,res,N,update);
    for(int i=0;i<N;i++) p[i]+=0.2*update[i];
  }while(iter++<maxiter);
  free(jac);
  free(res);
  free(update);
  free(Vnorm);
  free(Vprevnorm);
}
// J[i] is current density at i+0.5.
// Normalized V is passed. 
void compute_J(double *J, double *V, double *n, double u, int N) {
  for (int i = 0; i < N - 1; i++) {
    if(IN_OX(i)) //Inside Oxide Current is zero
      J[i]=0;
    else //In bulk of semiconductor, this the equation
      J[i] = (kB * mos.T *u / mos.dx)* (n[i+1]*B(V[i+1]-V[i])-n[i]*B(-(V[i+1]-V[i])));
  }
  J[N-1]= (kB * mos.T *u/mos.dx)*((n_teq)-n[N-1]); //Ohmic Boundary at endpoint
  return;
}

// J[i] is current density at i+0.5.
// Normalized V is passed. 
void compute_Jp(double *J, double *V, double *p, double u, int N) {
  for (int i = 0; i < N -1; i++) {
    if(IN_OX(i)) //Inside Oxide Current is zero
      J[i]=0;
    else //In bulk of semiconductor, this the equation
      J[i] = -(kB * mos.T *u / mos.dx)* (p[i+1]*B(-(V[i+1]-V[i]))-p[i]*B(V[i+1]-V[i]));
  }
  J[N-1]= -(kB * mos.T *u/mos.dx)*((p_teq)-p[N-1]); //Ohmic Boundary at endpoint
  //In mine its not negated. in namans code it is negated.
  return;
}

// function that takes the matrices needed and throws back the jacobi of it
// rn this one is only for N but can be simply copied for P 
// Normalized V is passed
// u is the mobility
void computeJacobi_n(double *Jac,double u,double *V,double *p, int N){
  // so fill the jacobi with the required derivatives
  // jacobi is Nx3 as tridiagonal and the terms will come using the other matrices
  double T=mos.T;
  double dx=mos.dx;
  double dt=sim.dt;
  double C = (kB * T * u) / (2.0 * q * dx * dx);
  for(int i=0;i<N;i++){

    if(IN_OX(i)||i==0){ //inside oxide conc is unchanging. Thus jacobian is such. Residual is zero respectively
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;
    }
    else if(IN_OX(i-1)){ //At Oxide-Semiconductor boundary, J inside is zero. so some terms are not there
      Jac[3*i]=0;
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*p[i]/2 + C*B(V[i]-V[i+1])+1/dt;
      // DRi/Dni+1
      Jac[3*i+2]=-C*B(-V[i]+V[i+1]) ;
    }
    else if(i==N-1){ //At endpoint(Ohmic boundary condition)
      /*Jac[3*i]=-C*B(V[i-1]-V[i]);
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*p[i]/2 + C*B(-V[i-1]+V[i]) + C+1/dt; //limiting case where the thing becomes 1 for constant potential
      // DRi/Dni+1
      Jac[3*i+2]=0;*/
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;

    }
    else{ //In bulk
      Jac[3*i]=-C*B(V[i-1]-V[i]);
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*p[i]/2 + C*B(-V[i-1]+V[i]) + C*B(V[i]-V[i+1])+1/dt;
      // DRi/Dni+1
      Jac[3*i+2]=-C*B(-V[i]+V[i+1]) ;
    }

  }
}

void computeJacobi_p(double *Jac,double u,double *V,double *n, int N){
  // so fill the jacobi with the required derivatives
  // jacobi is Nx3 as tridiagonal and the terms will come using the other matrices
  double T=mos.T;
  double dx=mos.dx;
  double dt=sim.dt;
  double C = (kB * T * u) / (2.0 * q * dx * dx);
  for(int i=0;i<N;i++){

    if(IN_OX(i)||i==0){ //inside oxide conc is unchanging. Thus jacobian is such. Residual is zero respectively
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;
    }
    else if(IN_OX(i-1)){ //At Oxide-Semiconductor boundary, J inside is zero. so some terms are not there
      Jac[3*i]=0;
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*n[i]/2 - C*B(-V[i]+V[i+1])+1/dt;
      // DRi/Dni+1
      Jac[3*i+2]=C*B(V[i]-V[i+1]) ;
    }
    else if(i==N-1){ //At endpoint(Ohmic boundary condition)
      /*Jac[3*i]=C*B(-V[i-1]+V[i]);
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*n[i]/2 - C*B(V[i-1]-V[i]) - C+1/dt; //limiting case where the thing becomes 1 for constant potential
      // DRi/Dni+1
      Jac[3*i+2]=0;
      */
      Jac[3*i]=0;
      Jac[3*i+1]=1;
      Jac[3*i+2]=0;
    }
    else{ //In bulk
      Jac[3*i]=-C*B(-V[i-1]+V[i]);
      // DRi/Dni
      Jac[3*i+1]=mos.C_Rr*n[i]/2 + C*B(V[i-1]-V[i]) + C*B(-V[i]+V[i+1])+1/dt;
      // DRi/Dni+1
      Jac[3*i+2]=-C*B(V[i]-V[i+1]) ;
    }

  }
}

void residual_n(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N){
  // function to return the values of residual
  // used crank nicolson scheme for this (can be changed if not working)
  // we need current density for this
  double *J = malloc(N*sizeof(double));
  compute_J(J,V,n,u,N);
  double *Jprev = malloc(N*sizeof(double));
  compute_J(Jprev,Vprev,nprev,u,N);

  for(int i=0;i<N;i++){
    // inside oxide and metal assumed no change in n and p and took values as zero
    if(IN_OX(i)||i==0){
      res[i]=0;
    }
    else if(i==N-1){
      res[i]=n[i]-n_teq;
    }
    else{
      // sign of electron is handled
      // WARNING: the scheme changes for the function definition of compute_J
      res[i]=(n[i]-nprev[i])/(sim.dt) - ((J[i]-J[i-1]+Jprev[i]-Jprev[i-1])/(2*q*mos.dx)) -(mos.Gr - mos.C_Rr*(n[i]*p[i]/2+nprev[i]*pprev[i]/2-mos.ni*mos.ni));
    }
  }
  free(J);
  free(Jprev);
}

void residual_p(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N){
  // function to return the values of residual
  // used crank nicolson scheme for this (can be changed if not working)
  // we need current density for this
  double *J = malloc(N*sizeof(double));
  compute_Jp(J,V,p,u,N);
  double *Jprev = malloc(N*sizeof(double));
  compute_Jp(Jprev,Vprev,pprev,u,N);

  for(int i=0;i<N;i++){
    // inside oxide and metal assumed no change in n and p and took values as zero
    if(IN_OX(i)||i==0){
      res[i]=0;
    }
    else if(i==N-1){
      res[i]=p[i]-p_teq;
    }
    else{
      // sign of electron is handled
      // WARNING: the scheme changes for the function definition of compute_J
      res[i]=(p[i]-pprev[i])/(sim.dt) + ((J[i]-J[i-1]+Jprev[i]-Jprev[i-1])/(2*q*mos.dx)) -(mos.Gr - mos.C_Rr*(n[i]*p[i]/2+nprev[i]*pprev[i]/2-mos.ni*mos.ni));
    }
  }
  free(J);
  free(Jprev);
}


