#include "parameter_fetch.h"
#include <stdio.h>

struct parameter mos;

void init_default_parameters(void) {
    mos.t_oxide = 5e-9;
    mos.area    = 1e-8;
    mos.height  = 1e-6;
    mos.nz=100;

    mos.eps_oxide = 3.9 * 8.854e-12;
    mos.eps_si    = 11.7 * 8.854e-12;
    mos.Na = 1e23;
    mos.Nd = 0;
    mos.mu_n = 0.135;
    mos.mu_p = 0.045;

    mos.Vg  = 0.0;
    mos.Vfb = -0.2;
    mos.Vth = 0.7;
    mos.T   = 300.0;
}
