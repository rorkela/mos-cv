
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Forward declaration
double* poisson(double *permittivity, double *conc,
                double Vapp1, double Vapp2, int n);

int main(void)
{
    int n = 10;              // grid size
    double Vapp1 = 0.0;      // left boundary voltage
    double Vapp2 = 1.0;      // right boundary voltage

    // allocate arrays
    double *permittivity = malloc(n * sizeof(double));
    double *conc = malloc(n * sizeof(double));

    // fill with simple test data
    for (int i=0; i<n; i++) {
        permittivity[i] = 1.0;     // constant ε
        conc[i] = 0.0;             // no charge inside
    }

    // call solver
    double *V = poisson(permittivity, conc, Vapp1, Vapp2, n);

    // print results
    printf("Node\tPotential (V)\n");
    for (int i=0; i<n; i++) {
        printf("%2d\t%.6f\n", i, V[i]);
    }

    // cleanup
    free(permittivity);
    free(conc);
    free(V);

    return 0;
}
