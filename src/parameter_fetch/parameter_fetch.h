#ifndef PARAMETER
#define PARAMETER
#define kB 1.380649e-23 // J/K
#define q 1.602176634e-19
#define IN_OX(i) ((i) * mos.dx <= mos.t_oxide)
#define n_teq mos.ni *mos.ni / mos.Na
#define p_teq mos.Na
// To store the input parameters for voltage.
// V=V.bias+V.sin*sin(2*M_PI*f*T);
struct signal {
  double bias; // Voltage bias
  double sin;  // Amplitude of sine
  double f;    // Frequency of sine
};
struct mos_param {
  double t_oxide; // Thickness of oxide. Again not needed ig
  double area;    // Area of mos(not needed ig. professor said unit area)
  double t_semi;  // Height of mesh
  int nz;         // Meshing points.
  double dx;      // Meshing distance. Computed. DONT SET.

  double eps_oxide; // Oxide permittivity
  double eps_si;    // Semiconductor permittivity
  double Na;        // Acceptor Doping conc
  double Nd;        // Donor Doping conc
  double mu_n;      // Mobility of n in silicon
  double mu_p;      // Mobility of p in silicon
  double ni;        // Intrinsic Carrier Concentration
  double Nc;        // Effective DOS in CB
  double Nv;        // Effective DOS in VB

  double T;    // Temperature
  double Gr;   // Generation rate uniformly assumed
  double C_Rr; // Radiative Recombination formula

  // add more as needed...
};
struct sim_param {
  double *Nd;   // donor doping. To be computed
  double *Na;   // acceptor doping. To be computed
  double *perm; // permittivity. To be computed
  double dt;    // Timestep.
  double *x;    // x axis position of each meshing point. To be computed
  int tdiv;     // Deprecated
};
extern struct sim_param sim;
extern struct mos_param mos;
void init_default_parameters(void);
void load_parameters_from_file(const char *filename); // Doesnt input everything yet. Use default parameters
void print_parameters(void);
void init_params();
void save_parameters_to_file(const char *filename);
void load_or_create_parameters(const char *filename);
#endif
