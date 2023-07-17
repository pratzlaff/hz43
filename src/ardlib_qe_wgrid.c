#include <ardlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char* argv[])
{

  const char* detsubsys;
  const char* obsfile;
  double x, y;
  float wmin, wmax;
  unsigned long i, n;

  Ardlib_Det_QE_Type *q;
  float *e, *qe;
  double winc;

  if (argc != 8) {
    fprintf(stderr, "Usage: %s detsubsys obsfile x y wmin wmax n\n", argv[0]);
    exit(1);
  }

  detsubsys = argv[1];
  obsfile = argv[2];
  x = atof(argv[3]);
  y = atof(argv[4]);
  wmin = atof(argv[5]);
  wmax = atof(argv[6]);
  n = atol(argv[7]);

  if (n > 1e6) {
    fprintf(stderr, "n = %lu is too large, exiting\n", n);
    exit(1);
  }

  if (wmin >= wmax) {
    fprintf(stderr, "wmin=%f >= wmax=%f, exiting\n", wmin, wmax);
    exit(1);
  }

  if (-1 == ardlib_initialize ("CHANDRA", obsfile))
     {
        fprintf (stderr, "Failed to initialize ardlib\n");
        exit (1);
     }

  /* generate energies */
  e = malloc(n * sizeof(float));
  qe = malloc(n * sizeof(float));
  if (!e || !qe) {
    fprintf(stderr, "could not allocate memory, exiting\n");
    exit(1);
  }

  winc = (wmax - wmin) / (n-1);
  for (i=0; i<n; ++i)
    e[i] = 12.39854/(wmax - i*winc);

  q = ardlib_open_det_qe(detsubsys);
  if (!q) {
    fprintf(stderr, "ardlib_open_det_qe(%s) failed, exiting\n", detsubsys);
    exit(1);
  }

  if (ardlib_compute_det_qe(q, x, y, e, n, qe)) {
    fprintf(stderr, "ardlib_compute_det_qe() failed, exiting\n");
    exit(1);
  }

  for (i=0; i<n; ++i)
    printf("%.2f\t%.6g\n", 12.39854/e[i], qe[i]);

  if (ardlib_close_det_qe(q)) {
    fprintf(stderr, "ardlib_close_det_qe() failed, exiting\n");
    exit(1);
  }

  free(e);
  free(qe);

  return 0;
}
