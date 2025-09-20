#include "../parameter_fetch/parameter_fetch.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Declare q and kB if not already in your headers

// Forward declaration of your poisson function
void poisson(double *V, double *n, double *p, double *permittivity, double Vbound1, double Vbound2);

int main() {
  // --- initialize MOS parameters ---
  mos.nz = 21;   // number of grid points
  mos.dx = 1e-7; // grid spacing 100 nm
  mos.T = 300;   // temperature
  mos.Nd = 1e23; // donor density
  mos.Na = 0;    // acceptor density (p-type doped MOSCAP example)

  int N = mos.nz;

  // --- allocate arrays ---
  double *V = calloc(N, sizeof(double));
  double *n = calloc(N, sizeof(double));
  double *p = calloc(N, sizeof(double));
  double *eps = calloc(N, sizeof(double));

  if (!V || !n || !p || !eps) {
    printf("Memory allocation failed!\n");
    return 1;
  }

  // --- initialize arrays ---
  for (int i = 0; i < N; i++) {
    n[i] = 0.0;               // no free electrons initially
    p[i] = 0.0;               // no free holes initially
    eps[i] = 11.7 * 8.85e-12; // silicon permittivity
    V[i] = 0.0;               // initial guess for potential
  }

  // --- boundary conditions ---
  double Vleft = 0.0;
  double Vright = 0.0;

  // --- call Poisson solver ---
  poisson(V, n, p, eps, Vleft, Vright);

  // --- print solution ---
  printf("x (m) \t V (V)\n");
  for (int i = 0; i < N; i++) {
    printf("%e \t %e\n", i * mos.dx, V[i]);
  }

  // --- free memory ---
  free(V);
  free(n);
  free(p);
  free(eps);

  return 0;
}
