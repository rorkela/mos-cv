#include "thomas.h"
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
