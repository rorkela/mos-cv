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
