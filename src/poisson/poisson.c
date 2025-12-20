#include "poisson.h"
#define MAX_ITER 20
/* This one will use newton rhapson to solve.
 * Takes input as
 * V - Voltage array and with initial guess supplied.
 * n - conc at that time instant
 * p - conc at that time instant
 * permittivity - permittivity of medium
 * Vbound1 and Vbound2 - Drichlet boundary
 * Updates V till convergence. After convergence exits.*/
void poisson(double *V, double *n, double *p, double Vbound1, double Vbound2) {
  int N = mos.nz;
  double *permittivity = sim.perm;
  double *Nd = sim.Nd;
  double *Na = sim.Na;
  double dx = mos.dx;
  int i;
  int iter = 0;
  // NOTE: Think of Jx=F where J is jacobian, X is update, F is -F(x) (- is
  // embedded inside beforehand for convenience) Since J is tridiagonal, i store
  // it in a Nx3 matrix. Saves time and space.
  double *J = malloc(3 * N * sizeof(double));
  double *X = malloc(N * sizeof(double));
  double *F = malloc(N * sizeof(double));
  // Jacobian
  // NOTE: Needs to be built once. For one run, the Jacobian is constant.
  J[0 * 3 + 1] = 1;
  J[(N - 1) * 3 + 1] = 1;
  J[0 * 3 + 0] = J[0 * 3 + 2] = 0;
  J[(N - 1) * 3 + 0] = J[(N - 1) * 3 + 2] = 0;
  for (i = 1; i < N - 1; i++) {
    J[i * 3] = (permittivity[i] + permittivity[i - 1]) / (2 * dx * dx);
    J[i * 3 + 1] = -(2 * permittivity[i] + permittivity[i + 1] + permittivity[i - 1]) / (2 * dx * dx) -
                   q * q * (p[i] + n[i]) / (kB * mos.T);
    J[i * 3 + 2] = (permittivity[i] + permittivity[i + 1]) / (2 * dx * dx);
  }
  do {
    // Setting boundary conditions
    V[0] = Vbound1;
    V[N - 1] = Vbound2;
    // Residual
    F[0] = 0;
    F[N - 1] = 0;
    for (int i = 1; i < N - 1; i++) {
      F[i] = (1 / (2 * dx * dx) *
                  ((permittivity[i] + permittivity[i - 1]) * V[i - 1] -
                   (2 * permittivity[i] + permittivity[i + 1] + permittivity[i - 1]) * V[i] +
                   (permittivity[i] + permittivity[i + 1]) * V[i + 1]) +
              q * (Nd[i] - Na[i] - n[i] + p[i]));
    }
    // Inverting Residual
    for (int i = 0; i < N; i++)
      F[i] = -F[i];
    thomas(J, F, N, X);
    // Update
    //  NOTE: Not accounting for boundary conditions. So loop starts at 1 and ends at N-2
    for (int i = 1; i < N - 1; i++)
      V[i] += X[i];
  } while (iter++ <= MAX_ITER);
  free(J);
  free(X);
  free(F);
}
