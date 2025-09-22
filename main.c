#include "main.h"

int main() {
  init_default_parameters();
  init_params();
  int N = mos.nz;

  // Defining
  struct signal Vin;
  int i;
  int dcdiv=10;
  double Vstart=-2;
  double Vend=2;
  Vin.f=10000;
  Vin.sin=0.01;
  double *C = malloc(dcdiv*sizeof(double));
  double *bias =malloc(dcdiv*sizeof(double));
  for (int i=0;i<dcdiv;i++) {
    Vin.bias=i*Vstart+(dcdiv-i)*Vend;
    bias[i]=Vin.bias;
    C[i] = solve_c(Vin);
  }
  plotxy(bias,C,dcdiv);
  free(C);
  free(bias);
}
