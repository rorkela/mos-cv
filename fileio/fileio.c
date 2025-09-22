#include "../main.h"

void plotxy(double *x, double *y, int N) {
  FILE *gp = popen("gnuplot -persistent", "w");
  if (!gp) {
    perror("gnuplot");
    return;
  }

  // Tell gnuplot to expect inline data
  fprintf(gp, "plot '-' with linespoints title 'XY Data'\n");

  // Send x–y pairs
  for (int i = 0; i < N; i++) {
    fprintf(gp, "%e %e\n", x[i], y[i]);
  }
  fprintf(gp, "e\n"); // end of dataset

  fflush(gp); // push data to gnuplot
  pclose(gp); // close pipe
}
void printarr(double *x, int N) {
  for (int i = 0; i < N; i++)
    printf("%e\n", x[i]);
}

void plotstate(double *x, double *V, double *n, double *p)
{
    int N=mos.nz;
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        perror("gnuplot");
        return;
    }

    // Use multiplot: 3 rows, 1 column
    fprintf(gp, "set multiplot layout 3,1 title 'State Plots'\n");

    // First subplot: V(x)
    fprintf(gp, "plot '-' with linespoints title 'V(x)'\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], V[i]);
    }
    fprintf(gp, "e\n");

    // Second subplot: n(x)
    fprintf(gp, "plot '-' with linespoints title 'n(x)'\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], n[i]);
    }
    fprintf(gp, "e\n");

    // Third subplot: p(x)
    fprintf(gp, "plot '-' with linespoints title 'p(x)'\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], p[i]);
    }
    fprintf(gp, "e\n");

    // End multiplot
    fprintf(gp, "unset multiplot\n");

    fflush(gp);
    pclose(gp);
}
