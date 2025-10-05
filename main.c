#include "main.h"
struct array_memory mem;
int main() {
  init_default_parameters();
  init_params();
  init_memory();

  // Defining
  struct signal Vin;
  int i;
  int dcdiv=31;
  double Vstart=2;
  double Vend=-1;
  Vin.f=10;
  Vin.sin=0.01;
  double *C = malloc(dcdiv*sizeof(double));
  double *bias =malloc(dcdiv*sizeof(double));
  FILE *out=fopen("charge.csv","w");
  for (i=0;i<dcdiv;i++) {
    printf("------%d/%d-------\n",i,dcdiv);
    Vin.bias=i*Vend/dcdiv+(dcdiv-i)*Vstart/dcdiv;
    bias[i]=Vin.bias;
    C[i] = solve_c(Vin);
    fprintf(out,"%e\t%e\n", bias[i],C[i]);
  }
  plotxy(bias,C,dcdiv);
  fclose(out);
  free(C);
  free(bias);
  free_memory();
}



void init_memory(){
  for(int i=0;i<10;i++)
  {
    mem.arrayn[i]=malloc(mos.nz*sizeof(double));
  }
  mem.array3n=malloc(3*mos.nz*sizeof(double));
}
void free_memory(){
  for(int i=0;i<10;i++)
  {
    free(mem.arrayn[i]);
  }
  free(mem.array3n);
}
