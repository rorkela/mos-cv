#ifndef PARAMETER
#define PARAMETER
struct parameter {
    double t_oxide;
    double area;
    double height;
    int nz;

    double eps_oxide;
    double eps_si;
    double Na;
    double Nd;
    double mu_n;
    double mu_p;

    double Vg;
    double Vfb;
    double Vth;
    double T;

    // add more as needed...
};
extern struct parameter mos;
void init_default_parameters(void);
#endif
