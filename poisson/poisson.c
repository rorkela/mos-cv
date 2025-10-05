#include "../main.h"
#define MAX_ITER 10
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
  double *J = mem.array3n;
  double *X = mem.arrayn[0];
  double *F = mem.arrayn[1];
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
}
void thomas(double *A, double *B, int N, double *x) {
  // here assumed the A matrix passed to this is already compacted tridiagonal Nx3
  // B matrix should be Nx1
  // X is the update so Nx1

  // pass the number of rows in this function which is just mos.nz
  double *c_new = malloc(N * sizeof(double)); // the diagonal elements above the one wala in the T matrix
  double *d_new = malloc(N * sizeof(double)); // the B matrix after the row operations and normalisation
  // Forward elimination
  // normalise the first and put the values already
  c_new[0] = A[2] / A[1]; // c1 / b1
  d_new[0] = B[0] / A[1]; // d1 / b1

  for (int i = 1; i < N; i++) {
    double denominator = A[i * 3 + 1] - A[i * 3 + 0] * c_new[i - 1];
    c_new[i] = (i == N - 1) ? 0.0 : A[i * 3 + 2] / denominator;
    d_new[i] = (B[i] - A[i * 3 + 0] * d_new[i - 1]) / denominator;
  }
  // Back substitution
  x[N - 1] = d_new[N - 1];
  for (int i = N - 2; i >= 0; i--) {
    x[i] = d_new[i] - c_new[i] * x[i + 1];
  }
  free(c_new);
  free(d_new);
}
