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
void poisson(double *V, double *n, double *p, 
             double Vbound1, double Vbound2) {
  int N = mos.nz;
  double *permittivity=sim.perm;
  double *Nd=sim.Nd;
  double *Na=sim.Na;
  double dx = mos.dx;
  int i, j;
  int iter = MAX_ITER;
  // NOTE: Think of Jx=F where J is jacobian, X is update, F is -F(x) (- is
  // embedded inside beforehand for convenience) Since J is tridiagonal, i store
  // it in a Nx3 matrix. Saves time and space.
  double *J = calloc(N * 3, sizeof(double));
  double *X = malloc(N * sizeof(double));
  double *F = malloc(N * sizeof(double));
  double *C = malloc(N * sizeof(double));
  double *D = malloc(N * sizeof(double));
  double Fnorm = 0;
  double Fnorm_prev = 0;
  // Initializting residual
  F[0] = 0;
  F[N - 1] = 0;
  for (int i = 1; i < N - 1; i++)
    F[i] = -(1 / (2 * dx * dx) *
                 ((permittivity[i] + permittivity[i - 1]) * V[i - 1] -
                  (2 * permittivity[i] + permittivity[i + 1] +
                   permittivity[i - 1]) *
                      V[i] +
                  (permittivity[i] + permittivity[i + 1]) * V[i + 1]) -
             q * (Nd[i] - Na[i] - n[i] + p[i]));
  for (i = 0; i < N; i++)
    Fnorm = fmax(Fnorm, fabs(F[i]));
  do {
    // Setting boundary conditions
    V[0] = Vbound1;
    V[N - 1] = Vbound2;
    J[0 * 3 + 1] = 1;
    J[(N - 1) * 3 + 1] = 1;
    J[0 * 3 + 0] = J[0 * 3 + 2] = 0;
    J[(N - 1) * 3 + 0] = J[(N - 1) * 3 + 2] = 0;
    for (i = 1; i < N - 1; i++) {
      J[i * 3] = (permittivity[i] + permittivity[i - 1]) / (2 * dx * dx);
      J[i * 3 + 1] =
          -(2 * permittivity[i] + permittivity[i + 1] + permittivity[i - 1]) /
              (2 * dx * dx) -
          q * q * (p[i] + n[i]) / (kB * mos.T);
      J[i * 3 + 2] = (permittivity[i] + permittivity[i + 1]) / (2 * dx * dx);
    }

    // NOTE: To solve this matrix Ax=B, LU decomposition can be used. However
    // this is tridiagonal so Thomas algorithm gives O(n) time complexity.
    C[0] = J[2] / J[1];
    D[0] = F[0] / J[1];
    for (i = 1; i < N; i++) {
      C[i] = J[i * 3 + 2] / (J[i * 3 + 1] - J[i * 3] * C[i - 1]);
      D[i] =
          (F[i] - J[i * 3] * D[i - 1]) / (J[i * 3 + 1] - J[i * 3] * C[i - 1]);
    }
    X[N - 1] = D[N - 1];
    for (i = N - 2; i >= 0; i--) {
      X[i] = D[i] - C[i] * X[i + 1];
    }
    for (int i = 1; i < N - 1; i++) // Not accounting for boundary conditions.
                                    // So loop starts at 1 and ends at N-2
    {
      V[i] += X[i];
    }
    F[0] = 0;
    F[N - 1] = 0;
    Fnorm_prev = Fnorm;
    for (int i = 1; i < N - 1; i++) {
      F[i] = -(1 / (2 * dx * dx) *
                   ((permittivity[i] + permittivity[i - 1]) * V[i - 1] -
                    (2 * permittivity[i] + permittivity[i + 1] +
                     permittivity[i - 1]) *
                        V[i] +
                    (permittivity[i] + permittivity[i + 1]) * V[i + 1]) -
               q * (Nd[i] - Na[i] - n[i] + p[i]));
    }
    // Computing norm for stopping in convergence
    Fnorm = 0;
    for (i = 0; i < N; i++)
      Fnorm = fmax(Fnorm, fabs(F[i]));

    if (Fnorm / Fnorm_prev < 0.01)
      break;
  } while (iter--);
  free(J);
  free(X);
  free(F);
  free(C);
  free(D);
}
