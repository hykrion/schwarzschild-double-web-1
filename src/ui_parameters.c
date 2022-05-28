#include "header/ui_parameters.h"

#include "header/globals.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* N�mero de par�metros + 1 */
#define UI_PARAM_SIZE 7
/* Longitud m�xima de los nombres de los par�metros */
#define UI_PARAM_LEN 16

#define UI_PARAM_A 0
#define UI_PARAM_B 1
#define UI_PARAM_LL 2
#define UI_PARAM_RS 3
#define UI_PARAM_WMIN 4
#define UI_PARAM_WMAX 5
#define UI_PARAM_NW 6

double uiParams[UI_PARAM_SIZE];

/**
  @brief  Leer los diferentes par�metros que modificaresmos desde la GUI
          Por simplicidad estos par�metros est�n en modo texto en un fichero,
          parameters.txt, que tendr� un par�metro por cada l�nea.
  */
void readParams(double uiparams[])
{
  char str[UI_PARAM_SIZE][UI_PARAM_LEN];

  int i;
  for(i = 0; i < UI_PARAM_SIZE; i++)
    memset(str[i], '\0', UI_PARAM_LEN);

  /* Leer los par�metros de fichero */
  FILE *fp;
  fp = fopen("parameters.txt", "r");

  i = 0;
  while(!feof(fp))
  {
    fgets(str[i], UI_PARAM_LEN - 1, fp);
    i++;
  }

  fclose(fp);

  i = 0;
  while(i < UI_PARAM_SIZE)
  {
    uiparams[i] = atof(str[i]);
    i++;
  }
}

double getA()
{
  return uiParams[UI_PARAM_A];
}

double getB()
{
  return uiParams[UI_PARAM_B];
}

double getH()
{
  double a = getA();
  double b = getB();
  double n = NODES - 1;

  return (b - a) / n;
}

int getN()
{
  return NODES - 1;
}

double getLL()
{
  return uiParams[UI_PARAM_LL];
}

double getRS()
{
  return uiParams[UI_PARAM_RS];
}

double getWMin()
{
  return uiParams[UI_PARAM_WMIN];
}

double getWMax()
{
  return uiParams[UI_PARAM_WMAX];
}

double getHW()
{
  double wm = getWMin();
  double wM = getWMax();
  double nW = uiParams[UI_PARAM_NW];

  return (wM - wm) / nW;
}

int getNW()
{
  return uiParams[UI_PARAM_NW];
}
