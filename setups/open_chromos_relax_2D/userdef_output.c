#include "pluto.h"

double VelocityRewriteCoeff(double x_loop);

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{

  int i, j, k, nv;
  double mu;
  double *v;
  double ***T, ***divB, ***av;

  v = ARRAY_1D(NVAR, double);

  T = GetUserVar("T");
  divB = GetUserVar("divB");
  av = GetUserVar("av");

  DOM_LOOP(k,j,i){
    NVAR_LOOP(nv) {
      v[nv] = d->Vc[nv][k][j][i];
      }
    mu = MeanMolecularWeight(v);
    T[k][j][i] = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * mu;
	  divB[k][j][i]	= 0.5*(EXPAND(
      + (d->Vc[BX1][k][j][i+1] - d->Vc[BX1][k][j][i-1]) / grid->dx[IDIR][i],
      + (d->Vc[BX2][k][j+1][i] - d->Vc[BX2][k][j-1][i]) / grid->dx[JDIR][j],
      + (d->Vc[BX3][k+1][j][i] - d->Vc[BX2][k-1][j][i]) / grid->dx[KDIR][k]
      ));
    av[k][j][i] = VelocityRewriteCoeff(grid->x[JDIR][j]);
  }

}
/* ************************************************************* */
void ChangeOutputVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;

}





