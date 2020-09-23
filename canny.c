#include <stdio.h> 

#include <math.h>

#include <stdlib.h>

#define PICSIZE 256
#define MAXMASK 100

int pic[PICSIZE][PICSIZE];
double outpic1[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
double final[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE];
double convy[PICSIZE][PICSIZE];
double ival[PICSIZE][PICSIZE];
int cand[PICSIZE][PICSIZE];
int histogram[PICSIZE];

main(argc, argv)
int argc;
char ** argv; {
  int i, j, p, q, s, t, mr, centx, centy;
  double xmaskval, ymaskval, sumx, sumy, sig, maxival, minval, maxval, ZEROTOL, slope, cutOff, percent;
  FILE * fo1, * fo2, * fo3, * fp1, * fopen();
  char * foobar;
  char throwaway[80];

  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      cand[i][j] = 0;
      final[i][j] = 0;
    }
  }

  argc--;
  argv++;
  foobar = * argv;
  fp1 = fopen(foobar, "rb");

  argc--;
  argv++;
  foobar = * argv;
  percent = atof(foobar);

  argc--;
  argv++;
  foobar = * argv;
  sig = atof(foobar);

  /* argc--; argv++;
  foobar = *argv; */
  fo1 = fopen("magnitude.pgm", "wb");

  fprintf(fo1, "P5\n");
  fprintf(fo1, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo1, "255\n");

  fo2 = fopen("peaks.pgm", "wb");

  fprintf(fo2, "P5\n");
  fprintf(fo2, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo2, "255\n");

  fo3 = fopen("final.pgm", "wb");

  fprintf(fo3, "P5\n");
  fprintf(fo3, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo3, "255\n");
  // Fix artifacts
  fgets(throwaway, 80, fp1);
  fgets(throwaway, 80, fp1);
  fgets(throwaway, 80, fp1);
  if ( !( (throwaway[0]=='2') && (throwaway[1]=='5') && (throwaway[2]=='5')))
  fgets(throwaway, 80, fp1); 

  mr = (int)(sig * 3);
  centx = (MAXMASK / 2);
  centy = (MAXMASK / 2);

  for (i = 0; i < 256; i++) {
    for (j = 0; j < 256; j++) {
      pic[i][j] = getc(fp1);
    }
  }

  for (p = -mr; p <= mr; p++) {
    for (q = -mr; q <= mr; q++) {
      /* maskval = ((2-(((p*p)+(q*q))/(sig*sig)))*
      (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));*/
      xmaskval = -q * (exp(-1 * (((p * p) + (q * q)) / (2 * (sig * sig)))));
      ymaskval = -p * (exp(-1 * (((p * p) + (q * q)) / (2 * (sig * sig)))));

      maskx[p + centy][q + centx] = xmaskval;
      masky[p + centy][q + centx] = ymaskval;
    }
  }

  for (i = mr; i <= 255 - mr; i++) {
    for (j = mr; j <= 255 - mr; j++) {
      sumx = 0;
      sumy = 0;
      for (p = -mr; p <= mr; p++) {
        for (q = -mr; q <= mr; q++) {
          sumx += pic[i + p][j + q] * maskx[p + centy][q + centx];
          sumy += pic[i + p][j + q] * masky[p + centy][q + centx];
        }
      }
      outpic1[i][j] = sumx;
      outpic2[i][j] = sumy;
      convx[i][j] = sumx;
      convy[i][j] = sumy;
    }
  }

  maxival = 0;
  for (i = mr; i < 256 - mr; i++) {
    for (j = mr; j < 256 - mr; j++) {
      ival[i][j] = sqrt((double)((convx[i][j] * convx[i][j]) +
        (convy[i][j] * convy[i][j])));
      if (ival[i][j] > maxival)
        maxival = ival[i][j];

    }
  }

  for (i = mr; i < PICSIZE - mr; i++) {
    for (j = mr; j < PICSIZE - mr; j++) {
      if (convx[i][j] == 0.0) {
        //continue;
        convx[i][j] = 0.00001;
      }

      slope = (convy[i][j] / convx[i][j]);

      if (slope <= 0.4142 && slope > -0.4142) {
        if (ival[i][j] > ival[i][j - 1] && ival[i][j] > ival[i][j + 1]) {
          cand[i][j] = 255;
        }
      } else if (slope <= 2.4142 && slope > 0.4142) {
        if (ival[i][j] > ival[i - 1][j - 1] && ival[i][j] > ival[i + 1][j + 1]) {
          cand[i][j] = 255;
        }
      } else if (slope <= -0.4142 && slope > -2.4142) {
        if (ival[i][j] > ival[i + 1][j - 1] && ival[i][j] > ival[i - 1][j + 1]) {
          cand[i][j] = 255;
        }
      } else {
        if (ival[i][j] > ival[i - 1][j] && ival[i][j] > ival[i + 1][j]) {
          cand[i][j] = 255;
        }
      }
    }
  }

  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      ival[i][j] = (ival[i][j] / maxival) * 255;
      fprintf(fo1, "%c", (char)((int)(ival[i][j])));

    }
  }

  fclose(fo1);
  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      fprintf(fo2, "%c", (char)((int) cand[i][j]));

    }
  }

  fclose(fo2);

  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      (histogram[(int) ival[i][j]]) ++;
    }
  }

  double total = 0;
  double area = percent * 0.01 * (PICSIZE * PICSIZE);

  for (i = PICSIZE; i > 0; i--) {
    total += histogram[i];
    if (total > (int) area)
      break;
  }

  double hi = i;
  double lo = hi * 0.35;

  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      if (cand[i][j] == 255) {
        if (ival[i][j] > hi) {
          final[i][j] = 255;
          cand[i][j] = 0;
        } else if (ival[i][j] < lo) {
          cand[i][j] = 0;
          final[i][j] = 0;
        }
      }
    }
  }

  int moretodo = 1;
  while (moretodo) {
    moretodo = 0;
    for (i = 0; i < PICSIZE; i++) {
      for (j = 0; j < PICSIZE; j++) {
        if (cand[i][j] == 255) {
          for (p = -1; p <= 1; p++) {
            for (q = -1; q <= 1; q++) {
              if (final[i + p][j + q] == 255) {
                final[i][j] = 255;
                cand[i][j] = 0;
                moretodo = 1;
              }
            }
          }
        }
      }
    }
  }

  for (i = 0; i < PICSIZE; i++) {
    for (j = 0; j < PICSIZE; j++) {
      fprintf(fo3, "%c", (char)((int)(final[i][j])));
    }
  }

  return 0;
}
