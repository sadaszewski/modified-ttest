/* $ID: mod_ttest.c, last updated 2018-08-07, F.Osorio */

#include <math.h>

/* some definitions */
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))
#define MIN(a,b)    (((a)<(b)) ? (a) : (b))
#define EQUAL(a,b)  (((a)!=(b)) ? (0) : (1))
#define SQR(x)      ((x) * (x))
#define repeat      for(;;)
#define HYPOT(dx, dy) (hypot(dx, dy))

/* dims structure */
typedef struct DIMS_struct {
  int
    n,              /* number of observations */
    p,              /* number of variables */
    nclass;         /* total of classes */
} DIMS_struct, *DIMS;

/* find_interval() */
int
find_interval(double *y, int n, double x)
{
  for (int i = 0; i < n; i++) {
    if (x <= y[i])
      return i;
  }
  return (n - 1);
}

/* distance_max() */
double
distance_max(double *xpos, double *ypos, int n)
{ /* computes the maximum distance among locations */
  double dx, dy, val = 0.0;

  for (int j = 0; j < n; j++) {
    for (int i = j + 1; i < n; i++) {
      dx = (xpos[j] - xpos[i]);
      dy = (ypos[j] - ypos[i]);
      val = MAX(val, hypot(dx, dy));
    }
  }
  return val;
}

/* set_bounds() */
void
set_bounds(DIMS dims, double maxdist, int do_half, double *upper_bounds)
{
  double accum = 0.0, half = 0.5, length;

  if (do_half)
    maxdist *= half;
  length = maxdist / dims->nclass;
  for (int i = 0; i < dims->nclass; i++) {
    accum += length;
    upper_bounds[i] = accum;
  }
}

/* static functions.. */
void MoranI(double *, double *, DIMS, double *, double *, double *, double *, double *);
double estimated_ESS(double *, double *, DIMS, double *, double *);
void mod_ttest(double *, double *, DIMS, double *, double *, double *, double *, double *, double *, double *);
/* ..end declarations */

void online_covariance(double *x, double *y, int n, double *xbar, double *ybar, double *xvar, double *yvar) {
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_x2 = 0.0;
    double sum_y2 = 0.0;
    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
    }
    *xbar = sum_x / n;
    *ybar = sum_y / n;
    *xvar = (sum_x2 / n) - SQR(*xbar);
    *yvar = (sum_y2 / n) - SQR(*ybar);
}

void
MoranI(double *x, double *y, DIMS dims, double *xpos, double *ypos, double *upper_bounds,
  double *card, double *index)
{ /* Moran's I */
  int pos;
  double dx, dy, distance, sx, sy, xbar, xvar, ybar, yvar, wts;

  online_covariance(x, y, dims->n, &xbar, &ybar, &xvar, &yvar);

  for (int k = 0; k < dims->nclass; k++) {
    sx = sy = wts = 0.0;
    for (int j = 0; j < dims->n; j++) {
      for (int i = j + 1; i < dims->n; i++) {
        dx = (xpos[i] - xpos[j]);
        dy = (ypos[i] - ypos[j]);
        distance = HYPOT(dx, dy);
        pos = find_interval(upper_bounds, dims->nclass, distance);
        if (pos == k) {
          wts++;
          sx += (x[i] - xbar) * (x[j] - xbar);
          sy += (y[i] - ybar) * (y[j] - ybar);
        }
      }
    }
    index[k] = (sx / wts) / xvar;
    index[k + dims->nclass] = (sy / wts) / yvar;
    card[k] = wts;
  }
}

double
estimated_ESS(double *xpos, double *ypos, DIMS dims, double *upper_bounds, double *imoran)
{
  int pos;
  double corx, cory, dx, dy, distance, rxx, ryy, sxx, syy, sxy, trxx, tryy, trxy, ans;

  /* initialization of correlation matrices */
  sxx = syy = sxy = trxy = 0.0;
  for (int j = 0; j < dims->n; j++) {
    rxx = ryy = 0.0;
    for (int i = 0; i < dims->n; i++) {
      if (i != j) {
        dx = (xpos[i] - xpos[j]);
        dy = (ypos[i] - ypos[j]);
        distance = HYPOT(dx, dy);
        pos = find_interval(upper_bounds, dims->nclass, distance);
        corx = imoran[pos];
        cory = imoran[pos + dims->nclass];
      } else
        corx = cory = 1.0;
      rxx  += corx;
      ryy  += cory;
      trxy += corx * cory;
    }
    sxx += rxx;
    syy += ryy;
    sxy += rxx * ryy;
  }

  /* computation of the effective sample size */
  trxx = (double) dims->n - sxx / dims->n;
  tryy = (double) dims->n - syy / dims->n;
  trxy += (sxx * syy / dims->n - 2.0 * sxy) / dims->n;
  ans = trxx * tryy / trxy;

  return (ans + 1.0);
}

void
mod_ttest(double *x, double *y, DIMS dims, double *xpos, double *ypos, double *upper_bounds,
  double *cor, double *card, double *imoran, double *stats)
{
  int lower_tail = 0, log_p = 0;
  double ESS, df, F, R; //, pval;

  MoranI(x, y, dims, xpos, ypos, upper_bounds, card, imoran);

  ESS  = estimated_ESS(xpos, ypos, dims, upper_bounds, imoran);
  R    = *cor; /* only for p = 2! */
  df   = ESS - 2.0;
  F    = SQR(R) / (1.0 - SQR(R)); /* unscaled F */
  // pval = pf(df * F, 1.0, df, lower_tail, log_p);

  /* save results */
  stats[0] = ESS;
  stats[1] = F;
  stats[2] = df;
}
