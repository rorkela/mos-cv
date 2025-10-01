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
    int N = mos.nz;
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        perror("gnuplot");
        return;
    }

    // Define datablocks first
    fprintf(gp, "$VDATA << EOD\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], V[i]);
    }
    fprintf(gp, "EOD\n");

    fprintf(gp, "$NDATA << EOD\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], n[i]);
    }
    fprintf(gp, "EOD\n");

    fprintf(gp, "$PDATA << EOD\n");
    for (int i = 0; i < N; i++) {
        fprintf(gp, "%e %e\n", x[i], p[i]);
    }
    fprintf(gp, "EOD\n");

    // Now use multiplot with datablocks
    fprintf(gp, "set multiplot layout 3,1 title 'State Plots'\n");
    fprintf(gp, "plot $VDATA with lines title 'V(x)'\n");
    fprintf(gp, "plot $NDATA with lines title 'n(x)'\n");
    fprintf(gp, "plot $PDATA with lines title 'p(x)'\n");
    fprintf(gp, "unset multiplot\n");

    fflush(gp);
    pclose(gp);
}
