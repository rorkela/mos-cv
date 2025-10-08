#ifndef MAINH
#define MAINH
#include "carrier_continuity/carrier.h"
#include "fileio/fileio.h"
#include "parameter_fetch/parameter_fetch.h"
#include "poisson/poisson.h"
#include "solve_c/solve_c.h"
#include "carrier_steady_state/carrier_ss.h"
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct array_memory {
	double *arrayn [10]; //index to select a array
	double *array3n; //Direct array for jacobian
};
extern struct array_memory mem;
void init_memory();
void free_memory();
#endif
