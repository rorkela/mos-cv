#include "main.h"
#include "parameter_fetch/parameter_fetch.h"
int main(int argc, char *argv[]) {
  // initialization
  init_default_parameters();
  init_params();
  if (argc < 3) {
    printf("Syntax: %s <parameter input file(if not found,it will be generated)> <output file>", argv[0]);
    return 0;
  }
  load_or_create_parameters(argv[1]);

  // Defining
  struct signal Vin;
  int i;
  int dcdiv = 31;
  double Vstart = 1.5;
  double Vend = -1;
  Vin.f = 1E9;
  // scanf("%lf",&mos.Gr);
  //  For output
  FILE *out = fopen(argv[2], "w");
  shit *output = malloc(dcdiv * sizeof(shit));

  // For progress printing IGNORE
  char *progress = malloc(dcdiv * sizeof(char) + 3);
  progress[0] = '[';
  progress[dcdiv + 1] = ']';
  progress[dcdiv + 2] = 0;
  for (int i = 0; i < dcdiv; i++)
    progress[i + 1] = ' ';
  // DONT IGNORE
  #pragma omp parallel for
  for (i = 0; i < dcdiv; i++) {
    printf("%s------%3d/%3d-------\n", progress, i, dcdiv);
    Vin.bias = i * Vend / dcdiv + (dcdiv - i) * Vstart / dcdiv;
    Vin.sin = (Vin.bias) / 10;
    output[i] = solve_c(Vin);
    progress[i + 1] = '=';
    printf("%s Bias=%e\tQdc=%e\n", progress, output[i].Vbias, output[i].Qdc);
  }
  // Numerical Derivative for DC Capacitance.
  for (i = 1; i < dcdiv; i++) {
    output[i].Cdc = (output[i].Qdc - output[i - 1].Qdc) / (output[i].Vbias - output[i - 1].Vbias);
  }
  output[0].Cdc = (output[1].Qdc - output[0].Qdc) / (output[1].Vbias - output[0].Vbias);
  for (int i = 0; i < dcdiv; i++)
    fprintf(out, "%e\t%e\t%e\t%e\t%e\t%e\n", output[i].Vbias, output[i].Qdc, output[i].dVac, output[i].dQac,
            output[i].Cac, output[i].Cdc);
  fclose(out);
  free(output);
}
