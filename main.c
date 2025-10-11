#include "main.h"
#include "parameter_fetch/parameter_fetch.h"
struct array_memory mem;
int main(int argc, char *argv[]) {
  //initialization
  init_default_parameters();
  init_params();
  if(argc<3)
  {
    printf("Syntax: %s <parameter input file(if not found,it will be generated)> <output file>", argv[0]);
    return 0;
  }
  load_or_create_parameters(argv[1]);


  // Defining
  struct signal Vin;
  int i;
  int dcdiv=31;
  double Vstart=1.5;
  double Vend=-1;
  Vin.f=1E9;
  //scanf("%lf",&mos.Gr);
  FILE *out=fopen(argv[2],"w");
  shit *output=malloc(dcdiv*sizeof(shit));
  #pragma omp parallel for
  for (i=0;i<dcdiv;i++) {
    printf("------%d/%d-------\n",i,dcdiv);
    Vin.bias=i*Vend/dcdiv+(dcdiv-i)*Vstart/dcdiv;
    Vin.sin=(Vin.bias)/10;
    output[i] = solve_c(Vin);
  }
  //Numerical Derivative for DC Capacitance.
  for (i=1;i<dcdiv;i++){
    output[i].Cdc=(output[i].Qdc-output[i-1].Qdc)/(output[i].Vbias-output[i-1].Vbias);
  }
  output[0].Cdc=(output[1].Qdc-output[0].Qdc)/(output[1].Vbias-output[0].Vbias);
  for(int i=0;i<dcdiv;i++)
    fprintf(out,"%e\t%e\t%e\t%e\t%e\t%e\n", output[i].Vbias,output[i].Qdc,output[i].dVac,output[i].dQac,output[i].Cac,output[i].Cdc);
  fclose(out);
  free(output);
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
