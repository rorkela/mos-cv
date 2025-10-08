#ifndef CARRIER_SS
#define CARRIER_SS
void carrier_continuity_ss(double *V, double *Vprev, double *nprev, double *pprev, double *n, double *p, int N);
void compute_J_ss(double *J, double *V, double *n, double u, int N);
void computeJacobi_n_ss(double *Jac,double u,double *V,double *p, int N);
void computeJacobi_p_ss(double *Jac,double u,double *V,double *n, int N);
void residual_n_ss(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N);
void residual_p_ss(double* res,double* n,double* p,double* nprev,double* pprev,double* V,double* Vprev,double u,int N);
#endif
