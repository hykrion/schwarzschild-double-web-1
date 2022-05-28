#include "header/txt_data.h"

#include "header/ui.h"
#include "header/phi.h"
#include "header/tortoise.h"

#include <stdio.h>
#include <math.h>

void
txt_data_wave(void)
{
  FILE *fp;
  fp = fopen("wave.txt", "w");

  double x = ui_get_a();
  double h = ui_get_h_forward();
  double* phiRea = phi_get_rea();
  double* phiImg = phi_get_img();
  double* phiV = phi_get_v();

  int i;
  for (i = 0; i < NODES; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", x, phiRea[i], phiImg[i], phiV[i]);
    x += h;
  }
  fclose(fp);
}

/**
  @brief  Compare the calculated value to the analytical one
*/
void
txt_data_turtle()
{
  FILE *fp;
  fp = fopen("turtle.txt", "w");

  /* Calculated value */
  TPoint2D *xy = tortoise_get_xy();
  /*
  Analytical value
  We're calculating y(x) but plotting x(y) so we must invert the function in the
  gnuplot file. See gui.tcl
  */
  tortoise_calculate_yx_analytical();
  TPoint2D *yxAnalytical = tortoise_get_yx_analytical();

  fprintf(fp, "#%s %34s %44s\n", "x", "y", "analitica");
  int i;
  for (i = 0; i < NODES; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\t\n", point2D_get_x(xy[i]), point2D_get_y(xy[i]), point2D_get_y(yxAnalytical[i]));
  }
  fclose(fp);
}

void
txt_data_rt()
{
  FILE *fp;
  fp = fopen("coefficients.txt", "w");

  int nW = ui_get_nW();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double *T = phi_get_T();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    double err = log10(fabs(1 - R[i] - T[i]));
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\t%.32f\n", w, R[i], T[i], R[i] + T[i], err);
  }
  fclose(fp);
}

void
txt_data_sigma()
{

  FILE *fp;
  fp = fopen("sigma-l.txt", "w");

  double wMin = ui_get_wMin();
  int nW = ui_get_nW();
  double hW = ui_get_hW();
  double *deltaL = phi_get_deltaL();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    fprintf(fp, "%.32f\t%.32f\n", w, deltaL[i]);
  }
  fclose(fp);
}
