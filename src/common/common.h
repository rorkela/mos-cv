#ifndef COMMON
#define COMMON
//Define all structs and constants here
#define IN_OX(i) ((i) * mos.dx <= mos.t_oxide)
#define n_teq mos.ni *mos.ni / mos.Na
#define p_teq mos.Na

extern const double kB; // j/k
extern const double q; 
typedef struct {
  double Qdc;   // Charge in DC
  double Vbias; // DC Bias
  double dQac;  // Charge for AC signal
  double dVac;  // Amplitude of the AC signal
  double Cdc;   // Capacitance for dc
  double Cac;   // Capacitance for ac
} outputarr;

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
extern struct mos_param mos;
extern struct sim_param sim;
#endif
