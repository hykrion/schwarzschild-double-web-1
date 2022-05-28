#ifndef TORTOISE_H_INCLUDED
#define TORTOISE_H_INCLUDED

#include "point2D.h"

struct tortoise_xyParams
{
  double rS;
};

/* Integration */
void tortoise_calculate_xy_forward(void);
void tortoise_calculate_xy_backguard(void);
int tortoise_xy_integration(double a,
                            double b,
                            int nodes,
                            double ci,
                            void *param);
/* x(y) */
TPoint2D* tortoise_get_xy(void);
double tortoise_get_xy_y(double x);
/* y(x) */
TPoint2D* tortoise_get_yx_analytical(void);
void tortoise_calculate_yx_analytical(void);
/* Utils */
void tortoise_reverse(TPoint2D* arr);

#endif // TORTOISE_H_INCLUDED
