#include "main.h"


int main()
{
    init_default_parameters();
    init_params();
    int N=mos.nz;

    // Defining 
    struct signal Vin;
    for(Vin.bias=-1;Vin.bias<=2;Vin.bias+=0.1)
    {
        double c=solve_c(Vin);
        printf("%e %e\n",Vin.bias,c);
    }

}

