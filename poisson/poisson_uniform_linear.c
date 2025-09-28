#include "../parameter_fetch/parameter_fetch.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define map(x, y) (x) * n + (y)
double *poisson(double *conc, double Vapp1, double Vapp2) {
  int n = mos.nz;
  double dx = mos.dx;
  int i, j;
  // NOTE: Think of matrix Ax=B where x=potentials
  double *A = calloc(n * 3, sizeof(double));
  A[0 * 3 + 1] = 1;
  A[(n - 1) * 3 + 1] = 1;
  for (i = 1; i < n - 1; i++) {
    A[i * 3] = (mos.eps_si) / (dx * dx);
    A[i * 3 + 1] = -(2 * mos.eps_si) / (dx * dx);
    A[i * 3 + 2] = (mos.eps_si) / (dx * dx);
  }
  double *B = malloc(n * sizeof(double));
  B[0] = Vapp1;
  B[n - 1] = Vapp2;
  for (int i = 1; i < n - 1; i++)
    B[i] = conc[i];
  // NOTE: To solve this matrix Ax=B, LU decomposition can be used. However this is tridiagonal so Thomas algorithm
  // gives O(n) time complexity.
  double *C = malloc(n * sizeof(double));
  double *D = malloc(n * sizeof(double));
  C[0] = A[2] / A[1];
  D[0] = B[0] / A[1];
  for (i = 1; i < n; i++) {
    C[i] = A[i * 3 + 2] / (A[i * 3 + 1] - A[i * 3] * C[i - 1]);
    D[i] = (B[i] - A[i * 3] * D[i - 1]) / (A[i * 3 + 1] - A[i * 3] * C[i - 1]);
  }
  free(B);
  free(A);
  double *X = malloc(n * sizeof(double));
  X[n - 1] = D[n - 1];
  for (i = n - 2; i >= 0; i--) {
    X[i] = D[i] - C[i] * X[i + 1];
  }
  return X;
}


