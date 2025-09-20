#ifndef SOLVE
#define SOLVE
double solve_c(struct signal Vin);      // To call other functions and find C for that signal
double solve_charge_density(double *V); // To solve for charge density
void copy_arr(double *source, double *target, int N);
void compute_delta(double *delta, double val, double valprev);
#endif
