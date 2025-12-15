#include "../main.h"
struct mos_param mos;
struct sim_param sim;
void init_default_parameters(void) {
  mos.t_oxide = 5e-7;
  mos.area = 10e-6;
  mos.t_semi = 10e-6;
  mos.nz = 200;
  mos.eps_oxide = 3.9 * 8.854e-12;
  mos.eps_si = 11.68 * 8.854e-12;
  mos.Na = 1e20;
  mos.Nd = 0;
  mos.ni = 1e16;
  mos.mu_n = 0.135;
  mos.mu_p = 0.045;
  mos.Nv = 1.8e25;
  mos.Nc = 3.2e25;

  mos.T = 300.0;
  mos.Gr = 1e30;
  mos.C_Rr = 1.1e-8;
}

void load_parameters_from_file(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (!f) {
    perror("Error opening parameter file");
    exit(1);
  }

  char key[64];
  double val;

  while (fscanf(f, "%63s %lf", key, &val) == 2) {
    if (strcmp(key, "t_oxide") == 0)
      mos.t_oxide = val;
    else if (strcmp(key, "area") == 0)
      mos.area = val;
    else if (strcmp(key, "t_semi") == 0)
      mos.t_semi = val;
    else if (strcmp(key, "nz") == 0)
      mos.nz = (int)val;
    else if (strcmp(key, "eps_oxide") == 0)
      mos.eps_oxide = val;
    else if (strcmp(key, "eps_si") == 0)
      mos.eps_si = val;
    else if (strcmp(key, "Na") == 0)
      mos.Na = val;
    else if (strcmp(key, "Nd") == 0)
      mos.Nd = val;
    else if (strcmp(key, "mu_n") == 0)
      mos.mu_n = val;
    else if (strcmp(key, "mu_p") == 0)
      mos.mu_p = val;
    else if (strcmp(key, "ni") == 0)
      mos.ni = val;
    else if (strcmp(key, "Nc") == 0)
      mos.Nc = val;
    else if (strcmp(key, "Nv") == 0)
      mos.Nv = val;
    else if (strcmp(key, "T") == 0)
      mos.T = val;
    else if (strcmp(key, "Gr") == 0)
      mos.Gr = val;
    else if (strcmp(key, "C_Rr") == 0)
      mos.C_Rr = val;
    else
      printf("Warning: Unknown parameter '%s' ignored\n", key);
  }

  fclose(f);
}

void save_parameters_to_file(const char *filename) {
  FILE *f = fopen(filename, "w");
  if (!f) {
    perror("Error creating parameter file");
    return;
  }

  fprintf(f, "t_oxide %e\n", mos.t_oxide);
  fprintf(f, "area %e\n", mos.area);
  fprintf(f, "t_semi %e\n", mos.t_semi);
  fprintf(f, "nz %d\n", mos.nz);
  fprintf(f, "eps_oxide %e\n", mos.eps_oxide);
  fprintf(f, "eps_si %e\n", mos.eps_si);
  fprintf(f, "Na %e\n", mos.Na);
  fprintf(f, "Nd %e\n", mos.Nd);
  fprintf(f, "mu_n %e\n", mos.mu_n);
  fprintf(f, "mu_p %e\n", mos.mu_p);
  fprintf(f, "ni %e\n", mos.ni);
  fprintf(f, "Nc %e\n", mos.Nc);
  fprintf(f, "Nv %e\n", mos.Nv);
  fprintf(f, "T %e\n", mos.T);
  fprintf(f, "Gr %e\n", mos.Gr);
  fprintf(f, "C_Rr %e\n", mos.C_Rr);

  fclose(f);
  printf("Default parameters saved to '%s'\n", filename);
}

void load_or_create_parameters(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (f) {
    fclose(f);
    load_parameters_from_file(filename);
    printf("Loaded parameters from '%s'\n", filename);
  } else {
    printf("File '%s' not found. Writing default parameters...\n", filename);
    save_parameters_to_file(filename);
  }
}

void init_params() {
  sim.perm = malloc(mos.nz * sizeof(double));
  sim.Na = malloc(mos.nz * sizeof(double));
  sim.Nd = malloc(mos.nz * sizeof(double));
  sim.x = malloc(mos.nz * sizeof(double));
  mos.dx = (mos.t_semi + mos.t_oxide) / (mos.nz - 1);
  sim.dt = 1;
  sim.tdiv = 512;
  for (int i = 0; i < mos.nz; i++) {
    if (IN_OX(i)) {
      sim.Na[i] = 0;
      sim.Nd[i] = 0;
      sim.perm[i] = mos.eps_oxide;
    } else {
      sim.Na[i] = mos.Na;
      sim.Nd[i] = mos.Nd;
      sim.perm[i] = mos.eps_si;
    }
    sim.x[i] = i * mos.dx;
  }
}
void print_parameters(void) {}
