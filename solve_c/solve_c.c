#include "../main.h"
#define MAX_ITER 40
double solve_c(struct signal Vin) {
  int N = mos.nz;
  double drichlet_factor = kB * mos.T * log(mos.Nc / mos.Nd); // WARNING: Hardcoded for n doped
  double delta = 0;
  int iter = 0;
  // Defining variables for charge density for C calculations
  double Qdc;   // For Q in DC steady state
  double Qac;   // For Maximum Q in AC
  double Qprev; // Temporary variable
  // Defining parameters for time
  int tstep = 0;
  int tstepmax = sim.tdiv * 5;
  sim.dt = 1 / Vin.f / sim.tdiv;
  // Initializing arrays for n p V and for previous time instant
  double *n = malloc(N * sizeof(double));        // n for present computations
  double *p = malloc(N * sizeof(double));        // p for present computations
  double *V = malloc(N * sizeof(double));        // V for present computations
  double *n_prev = malloc(N * sizeof(double));   // n in previous computation for same t
  double *p_prev = malloc(N * sizeof(double));   // p in previous computation for same t
  double *V_prev = malloc(N * sizeof(double));   // V in previous computation for same t
  double *n_prev_t = malloc(N * sizeof(double)); // n in previous time instant.
  double *p_prev_t = malloc(N * sizeof(double)); // p in previous time instant.
  double *V_prev_t = malloc(N * sizeof(double)); // V in previous time instant.
  
  // DEBUG:
  //FILE *chargedc=fopen("charge_trans_dc.txt","w");
  //FILE *chargeac=fopen("charge_trans_ac.txt","w");
  // Starting from thermal equilibrium conditions at t=0;
  for (int i = 0; i < N; i++) {
    if (IN_OX(i)) {
      n[i] = 0;
      p[i] = 0;
    }
    // WARNING: Assuming this is n doped. Change here if needed. Hardcoded for now
    else {
      n[i] = mos.Nd;
      p[i] = mos.ni * mos.ni / mos.Nd;
    }
    V[i] = Vin.bias * (1 - (double)i / (N - 1)) + drichlet_factor;
  }
  Qdc = solve_charge_density(V);
  while (tstep++ <= tstepmax) {
    iter=0;
    printf("iter=%d",iter);
    plotstate(sim.x,V,n,p);
    copy_arr(n, n_prev_t, N);
    copy_arr(p, p_prev_t, N);
    copy_arr(V, V_prev_t, N);
    do {

      copy_arr(n, n_prev, N);
      copy_arr(p, p_prev, N);
      copy_arr(V, V_prev, N);
      poisson(V, n, p, V[0], V[N - 1]);
      carrier_continuity(V, V_prev_t, n_prev_t, p_prev_t, n, p, N);
      // Logic for computing delta
      delta = 0;
      for (int i = 0; i < N; i++) {
        compute_delta(&delta, V[i], V_prev[i]);
        compute_delta(&delta, n[i], n_prev[i]);
        compute_delta(&delta, p[i], p_prev[i]);
      }
      // if (delta <= 1e-25)
      //   break;
    } while (iter++ <= MAX_ITER);
    Qprev = Qdc;
    Qdc = solve_charge_density(V);
    delta = fabs(1 - Qdc / Qprev);
    //fprintf(chargedc,"%e\n",Qdc);
    printf("%e\n",Qdc);
    //if (delta <= 5e-3)
    //  break; // Tolerance is 0.5% change
  }
  // TODO: plotstate(sim.x,V,n,p);
  printf("solve_c.c: Qdc=%e\n", Qdc);
  // AC Analysis
  tstep=0;
  tstepmax=sim.tdiv*5;
  Qac=Qdc;
  while (tstep++ <= tstepmax) {
    iter=0;
    copy_arr(n, n_prev_t, N);
    copy_arr(p, p_prev_t, N);
    copy_arr(V, V_prev_t, N);
    V[0]=Vin.bias+Vin.sin*sin(2*M_PI*(tstep/(double)sim.tdiv))+drichlet_factor;
    do {

      copy_arr(n, n_prev, N);
      copy_arr(p, p_prev, N);
      copy_arr(V, V_prev, N);
      poisson(V, n, p, V[0], V[N - 1]);
      carrier_continuity(V, V_prev_t, n_prev_t, p_prev_t, n, p, N);

      // Logic for computing delta
      delta = 0;
      for (int i = 0; i < N; i++) {
        compute_delta(&delta, V[i], V_prev[i]);
        compute_delta(&delta, n[i], n_prev[i]);
        compute_delta(&delta, p[i], p_prev[i]);
      }
      // if (delta <= 1e-25)
      //   break;
    } while (iter++ <= MAX_ITER);
    Qac = fmax(Qac,solve_charge_density(V));
    //fprintf(chargeac,"%e\n",solve_charge_density(V));
    printf("%d %d\n",iter,tstep);
  }
  double c=(Qac-Qdc)/Vin.sin;
  printf("Solved\n");
  //fclose(chargedc);
  //fclose(chargeac);
  free(n);
  free(p);
  free(V);
  free(n_prev);
  free(p_prev);
  free(V_prev);
  free(n_prev_t);
  free(p_prev_t);
  free(V_prev_t);
  return c;
}

double solve_charge_density(double *V) // To solve for charge density
{
  return -mos.eps_oxide * (V[2] - V[1]) / mos.dx;
}

void copy_arr(double *source, double *target, int N) {
  memcpy(target,source,N*sizeof(double));
}
void compute_delta(double *delta, double val, double valprev) { *delta = fmax(*delta, fabs(val - valprev)); }
