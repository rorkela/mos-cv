#ifndef CARRIER
#define CARRIER
void carrier_continuity(double *V, double *Vprev, double *nprev, double *pprev, double *n, double *p, int N);
void compute_J(double *J, double *V, double *n, double u, int N);
void computeJacobi_n(double *Jac,double u,double *V,double *p, int N);
void computeJacobi_p(double *Jac,double u,double *V,double *n, int N);
void residual_n(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N);
void residual_p(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N);
#endif
