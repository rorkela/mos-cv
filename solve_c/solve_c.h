#ifndef SOLVE
#define SOLVE
#include "../main.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../parameter_fetch/parameter_fetch.h"
double solve_c(struct signal Vin);// To call other functions and find C for that signal
void solve_charge_density(double *charge_density, double *n, double *p, int N); //To solve for charge density
#endif
