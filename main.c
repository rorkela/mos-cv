#include "main.h"

int main() {
  init_default_parameters();
  init_params();

  // Defining
  struct signal Vin;
  int i;
  int dcdiv=3;
  double Vstart=-1.9;
  double Vend=2;
  Vin.f=10;
  Vin.sin=0.01;
  double *C = malloc(dcdiv*sizeof(double));
  double *bias =malloc(dcdiv*sizeof(double));
  for (i=0;i<dcdiv;i++) {
    Vin.bias=i*Vend/dcdiv+(dcdiv-i)*Vstart/dcdiv;
    bias[i]=Vin.bias;
    C[i] = solve_c(Vin);
  }
  plotxy(bias,C,dcdiv);
  free(C);
  free(bias);
}
