#include "header/tortoise.h"

#include "header/ui.h"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#define TOR_SYS_DIM 1
#define TOR_MIN_INF -1000000

static TPoint2D m_xy[NODES];
static TPoint2D m_yxAnalytical[NODES];

/**
  @brief  xy Jacobian: rS/x(y)^2
*/
static int
jacoXY(double t,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  (void)(t);
  struct tortoise_xyParams *par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, rS/(y[0]*y[0]));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  IV ODE to solve x[y]' = 1 - rS/x(y)
*/
static int
funcXY(double t,
       const double y[],
       double dydt[],
       void* params)
{
  (void)(t);
  struct tortoise_xyParams* par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  dydt[0] = 1 - rS/y[0];

  return GSL_SUCCESS;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
/**
  @brief  Calculate coordinates 'xy' integrating forward
*/
void
tortoise_calculate_xy_forward(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double rS = ui_get_rS();

  struct tortoise_xyParams torParam = {rS};
  double ic = rS*(1 + exp(a/rS - 1));
  /* TODO You can experiment with other ic...*
  double delta = 1e-9;
  ic = rS + delta;
  ic = rS;
  */

  tortoise_xy_integration(a, b, NODES, ic, &torParam);
}

/**
  @brief  Calculate coordinates 'xy' integrating back-guard

  NOTE    We must leave the array as it was made integrating forward, from
          [-inf, inf]
*/
void
tortoise_calculate_xy_backguard(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double rS = ui_get_rS();

  struct tortoise_xyParams torParam = {rS};
  /*
  TODO
    why the IC is 'b - rS*log(b/rS - 1)' instead of 'b + rS*log(b/rS - 1)' that
    is the correct expression to y[x] or r*(r) ??
  */
  double ic = b - rS*log(b/rS - 1);
  /* TODO You can experiment with other ic...*
  ic = b;
  */
  tortoise_xy_integration(b, a, NODES, ic, &torParam);
  tortoise_reverse(m_xy);
}

int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ic,
                        void *param)
{
  int status = GSL_SUCCESS;

  struct tortoise_xyParams *par = (struct tortoise_xyParams*)(param);
  double rS = par->rS;
  int n = nodes - 1;
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct tortoise_xyParams params = {rS};

  gsl_odeiv2_system sys = {funcXY, jacoXY, TOR_SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ic};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] function value
    y[1] derivative function value
    */
    TPoint2D p2d = {x0, y[0]};
    m_xy[i] = p2d;
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

TPoint2D*
tortoise_get_xy(void)
{
  return m_xy;
}

/**
  @brief  Analytical yx result from [-inf, inf]
*/
void
tortoise_calculate_yx_analytical(void)
{
  double x = ui_get_a();
  double h = ui_get_h_forward();
  double rS = ui_get_rS();
  int i;
  for(i = 0; i < NODES; i++)
  {
    double yx = x + rS*log(x/rS - 1);

    // log(y <= 0) = -inf
    if(isnan(yx))
      yx = TOR_MIN_INF;

    TPoint2D p2d = {x, yx};
    m_yxAnalytical[i] = p2d;
    x += h;
  }
}

TPoint2D*
tortoise_get_yx_analytical(void)
{
  return m_yxAnalytical;
}

/**
  @brief  Get the 'y' part from 'xy'
*/
double
tortoise_get_xy_y(double x)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();

  /* We get a double but we need an int */
  int i = (x - a) / h;

  return point2D_get_y(m_xy[i]);
}

void
tortoise_reverse(TPoint2D* arr)
{
  TPoint2D tmp[NODES];

  int i;
  for(i = 0; i < NODES; i++)
    tmp[i] = arr[NODES - 1 - i];
  for(i = 0; i < NODES; i++)
    arr[i] = tmp[i];
}
