#include "main.h"

int main() {
  init_default_parameters();
  init_params();

  // Defining
  struct signal Vin;
  int i;
  int dcdiv=100;
  double Vstart=-1.9;
  double Vend=2;
  Vin.f=10;
  Vin.sin=0.01;
  double *C = malloc(dcdiv*sizeof(double));
  double *bias =malloc(dcdiv*sizeof(double));
  FILE *out=fopen("charge.csv","w");
  for (i=0;i<dcdiv;i++) {
    Vin.bias=i*Vend/dcdiv+(dcdiv-i)*Vstart/dcdiv;
    bias[i]=Vin.bias;
    C[i] = solve_c(Vin);
    fprintf(out,"%e\t%e\n", bias[i],C[i]);
  }
  plotxy(bias,C,dcdiv);
  fclose(out);
  free(C);
  free(bias);
}
