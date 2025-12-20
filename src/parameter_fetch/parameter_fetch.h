#ifndef PARAMETER
#define PARAMETER
#include "../common/common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
// To store the input parameters for voltage.
// V=V.bias+V.sin*sin(2*M_PI*f*T);
void init_default_parameters(void);
void load_parameters_from_file(const char *filename); // Doesnt input everything yet. Use default parameters
void print_parameters(void);
void init_params();
void save_parameters_to_file(const char *filename);
void load_or_create_parameters(const char *filename);
#endif
