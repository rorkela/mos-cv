#ifndef PARAMETER
#define PARAMETER
#define kB 8.617333262e-5
#define q 1.6e-19
struct parameter {
    double t_oxide;//Thickness of oxide. Again not needed ig
    double area;//Area of mos(not needed ig. professor said unit area)
    double t_semi; //Height of mesh
    int nz; //Meshing points.
    double dx;

    double eps_oxide;
    double eps_si;
    double Na; //Doping conc
    double Nd; //Doping conc
    double mu_n;// Mobility of n in silicon
    double mu_p; //Mobility of p in silicon
    double ni;

    double Vg;
    double Vfb;
    double Vth;
    double T;//Temperature
    double Gr; //Generation rate uniformly assumed
    double C_Rr; //Radiative Recombination formula

    // add more as needed...
};
struct sim_arrays  {
    double * Nd; //donor doping
    double * Na; //acceptor doping
    double *perm; //permittivity
};
extern struct sim_arrays sim;
extern struct parameter mos;
void init_default_parameters(void);
int load_parameters_from_file(const char *fname); //Doesnt input everything yet. Use default parameters
void print_parameters(void); 
void init_params();
#endif
