#ifndef CARRIER
#define CARRIER
void compute_J(double *J, double *V, double *n,double u, int N);
void compute_F(double *F, double * V, double * Vprev, double *nprev, double *pprev, double * n, double * p, int N, int is_n);
#endif
