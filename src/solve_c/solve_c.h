#ifndef SOLVE
#define SOLVE
#include <math.h>
#include <stdio.h>
#include "../carrier_continuity/carrier.h"
#include "../poisson/poisson.h"
#include "../fileio/fileio.h"
#include "../parameter_fetch/parameter_fetch.h"
#include "../common/common.h"
#include <stdlib.h>
outputarr solve_c(double Vin);        // To call other functions and find C for that signal
double solve_charge_density(double *V); // To solve for charge density
void copy_arr(double *source, double *target, int N);
void compute_delta(double *delta, double val, double valprev);
#endif
