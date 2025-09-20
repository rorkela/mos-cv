#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define map(x, y) (x) * n + (y)
// This one assumes that permittivity is not constant and also some weird parameters are there.
// WARNING: NOT USED. Refer to poisson.c
double *poisson(double *permittivity, double *conc, double Vapp1, double Vapp2, int n, double dx) {
  int i, j;
  // NOTE: Think of matrix Ax=B where x=potentials
  double *A = calloc(n * n, sizeof(double));
  A[map(0, 0)] = 1;
  A[map(n - 1, n - 1)] = 1;
  for (i = 1; i < n - 1; i++) {
    A[map(i, i - 1)] = (permittivity[i] + permittivity[i - 1]) / (2 * dx * dx);
    A[map(i, i)] = -(2 * permittivity[i] + permittivity[i + 1] + permittivity[i - 1]) / (2 * dx * dx);
    A[map(i, i + 1)] = (permittivity[i] + permittivity[i + 1]) / (2 * dx * dx);
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
  C[0] = A[map(0, 1)] / A[map(0, 0)];
  D[0] = B[0] / A[map(0, 0)];
  for (i = 1; i < n; i++) {
    C[i] = A[map(i, i + 1)] / (A[map(i, i)] - A[map(i, i - 1)] * C[i - 1]);
    D[i] = (B[i] - A[map(i, i - 1)] * D[i - 1]) / (A[map(i, i)] - A[map(i, i - 1)] * C[i - 1]);
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
