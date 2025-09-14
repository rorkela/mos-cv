#include "parameter_fetch.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
struct parameter mos;
struct sim_arrays sim;
void init_default_parameters(void) {
    mos.t_oxide = 5e-9 ;
    mos.area    = 1e-8;
    mos.t_semi  = 1e-8;
    mos.nz=100;
    mos.eps_oxide = 3.9 * 8.854e-12;
    mos.eps_si    = 11.68 * 8.854e-12;
    mos.Na = 1e23;
    mos.Nd = 0;
    mos.ni=1.45e10;
    mos.mu_n = 0.135;
    mos.mu_p = 0.045;

    mos.Vg  = 0.0;
    mos.Vfb = -0.2;
    mos.Vth = 0.7;
    mos.T   = 300.0;
    mos.Gr = 0;
    mos.C_Rr=1.1e-14;
}

int load_parameters_from_file(const char *fname) {
    FILE *fp = fopen(fname, "r");
    if (!fp) {
        perror("Error opening parameter file");
        return -1;
    }

    char name[64];
    double value;

    while (fscanf(fp, "%63s %lf", name, &value) == 2) {
        if (strcmp(name, "t_oxide") == 0) mos.t_oxide = value;
        else if (strcmp(name, "area") == 0) mos.area = value;
        else if (strcmp(name, "height") == 0) mos.t_semi = value;
        else if (strcmp(name, "eps_oxide") == 0) mos.eps_oxide = value;
        else if (strcmp(name, "eps_si") == 0) mos.eps_si = value;
        else if (strcmp(name, "Na") == 0) mos.Na = value;
        else if (strcmp(name, "Nd") == 0) mos.Nd = value;
        else if (strcmp(name, "mu_n") == 0) mos.mu_n = value;
        else if (strcmp(name, "mu_p") == 0) mos.mu_p = value;
        else if (strcmp(name, "Vg") == 0) mos.Vg = value;
        else if (strcmp(name, "Vfb") == 0) mos.Vfb = value;
        else if (strcmp(name, "Vth") == 0) mos.Vth = value;
        else if (strcmp(name, "T") == 0) mos.T = value;
        else {
            printf("Warning: Unknown parameter '%s' ignored.\n", name);
        }
    }

    fclose(fp);
    return 0;
}
void init_params()
{
    sim.perm=malloc(mos.nz*sizeof(double));
    sim.Na=malloc(mos.nz*sizeof(double));
    sim.Nd=malloc(mos.nz*sizeof(double));
    mos.dx=(mos.t_semi+mos.t_oxide)/(mos.nz-1);
    for(int i=0;i<mos.nz;i++)
    {
        if(mos.dx*i<=mos.t_oxide) 
        {
            sim.Na[i]=0;
            sim.Nd[i]=0;
            sim.perm[i]=mos.eps_oxide;
        }
        else
        {
            sim.Na[i]=mos.Na;
            sim.Nd[i]=mos.Nd;
            sim.perm[i]=mos.eps_si;
        }
    }
}
void print_parameters(void) {
    printf("MOSCAP parameters:\n");
    printf("  t_oxide   = %g m\n", mos.t_oxide);
    printf("  area      = %g m^2\n", mos.area);
    printf("  height    = %g m\n", mos.t_semi);
    printf("  eps_oxide = %g F/m\n", mos.eps_oxide);
    printf("  eps_si    = %g F/m\n", mos.eps_si);
    printf("  Na        = %g m^-3\n", mos.Na);
    printf("  Nd        = %g m^-3\n", mos.Nd);
    printf("  mu_n      = %g m^2/Vs\n", mos.mu_n);
    printf("  mu_p      = %g m^2/Vs\n", mos.mu_p);
    printf("  Vg        = %g V\n", mos.Vg);
    printf("  Vfb       = %g V\n", mos.Vfb);
    printf("  Vth       = %g V\n", mos.Vth);
    printf("  T         = %g K\n", mos.T);
}
