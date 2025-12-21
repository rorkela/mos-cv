#include "solve_c.h"
#define MAX_ITER 40
#define MIN_TSTEP 50
outputarr solve_c(double Vin) {
  outputarr output;
  int N = mos.nz;

  double drichlet_factor = kB * mos.T / q * log(mos.Nc / n_teq); // WARNING: Hardcoded for n doped
  double delta = 0;
  int iter = 0;

  // Defining variables for charge density for C calculations
  double Qdc;       // For Q in DC steady state
  double Qprev = 0; // Temporary variable

  // Defining parameters for time
  int tstep = 0;
  int tstepmax = 5000;
  sim.dt = 2e-12;

  // Initializing arrays for n p V and for previous time instant
  double *n = malloc(N * sizeof(double));        // n for present computations
  double *p = malloc(N * sizeof(double));        // p for present computations
  double *V = malloc(N * sizeof(double));        // V for present computations
  double *n_prev_t = malloc(N * sizeof(double)); // n in previous time instant.
  double *p_prev_t = malloc(N * sizeof(double)); // p in previous time instant.
  double *V_prev_t = malloc(N * sizeof(double)); // V in previous time instant.

  // Adding dirichlet factor; taking into account fermi-level reference
  Vin= Vin+ drichlet_factor;

  for (int i = 0; i < N; i++) {
    if (IN_OX(i)) {
      n[i] = 0;
      p[i] = 0;
    } else {
      n[i] = n_teq;
      p[i] = p_teq;
    }
    V[i] = Vin- (double)i / (N - 1) * (Vin- drichlet_factor);
  }
  V[N - 1] = drichlet_factor;

  // **** DC ANALYSIS STARTS HERE ****
  Qdc = solve_charge_density(V);
  while (tstep++ <= tstepmax) {
    iter = 0;

    copy_arr(n, n_prev_t, N);
    copy_arr(p, p_prev_t, N);
    copy_arr(V, V_prev_t, N);

    do {
      poisson(V, n, p, V[0], V[N - 1]);
      carrier_continuity(V, V_prev_t, n_prev_t, p_prev_t, n, p, N);
    } while (iter++ <= MAX_ITER);

    Qprev = Qdc;
    Qdc = solve_charge_density(V);
    delta = fabs(1 - Qdc / Qprev);

    if (delta <= 2e-6 && tstep > MIN_TSTEP)
      break;
  }

  output.Qdc = Qdc;
  output.Vbias = Vin;


  // UNCOMMENT THIS FOR V n p PLOTS
  // plotstate(sim.x, V, n, p);

  free(n);
  free(p);
  free(V);
  free(n_prev_t);
  free(p_prev_t);
  free(V_prev_t);
  return (output);
}

// Function to solve for charge density using gauss's law
double solve_charge_density(double *V) { return -mos.eps_oxide * (V[2] - V[1]) / mos.dx; }

// Function to copy an array into an other one
void copy_arr(double *source, double *target, int N) { memcpy(target, source, N * sizeof(double)); }

// Function to compute absolute error b/w consecutive iterations
void compute_delta(double *delta, double val, double valprev) { *delta = fmax(*delta, fabs(val - valprev)); }
