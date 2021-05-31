#include "pluto.h"

#define CONST_gsun  (27.4e3 / (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH))

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

  /*
  int i, j, k, nv;
  double *v;
  double ***Jx1, ***Jx2, ***Jx3;
  double ***Sx1, ***Sx2, ***Sx3;
  double ***Fx1, ***Fx2, ***Fx3;
  #if RESISTIVITY != NO
  double *J, *eta;
  #endif
  double E;
  double L = g_inputParam[HALF_LOOP_L];

  Jx1 = GetUserVar("Jx1");
  Jx2 = GetUserVar("Jx2");
  Jx3 = GetUserVar("Jx3");

  Sx1 = GetUserVar("Sx1");
  Sx2 = GetUserVar("Sx2");
  Sx3 = GetUserVar("Sx3");

  Fx1 = GetUserVar("Fx1");
  Fx2 = GetUserVar("Fx2");
  Fx3 = GetUserVar("Fx3");

  DOM_LOOP(k,j,i) {

    NVAR_LOOP(nv) {
      v[nv] = d->Vc[nv][k][j][i];
      }

    // Current
    Jx1[k][j][i] = 0.5*(
      + (d->Vc[BX3][k][j+1][i] - d->Vc[BX3][k][j-1][i]) / grid->dx[JDIR][j]
      - (d->Vc[BX2][k+1][j][i] - d->Vc[BX2][k-1][j][i]) / grid->dx[KDIR][k]);
    Jx2[k][j][i] = 0.5*(
      - (d->Vc[BX3][k][j][i+1] - d->Vc[BX3][k][j][i-1]) / grid->dx[IDIR][i]
      + (d->Vc[BX1][k+1][j][i] - d->Vc[BX1][k-1][j][i]) / grid->dx[KDIR][k]);
    Jx3[k][j][i] = 0.5*(
      + (d->Vc[BX2][k][j][i+1] - d->Vc[BX2][k][j][i-1]) / grid->dx[IDIR][i]
      - (d->Vc[BX1][k][j+1][i] - d->Vc[BX1][k][j-1][i]) / grid->dx[JDIR][j]);

    // Resistivity
    #if RESISTIVITY != NO
    J[IDIR] = Jx1[k][j][i];
    J[JDIR] = Jx2[k][j][i];
    J[KDIR] = Jx3[k][j][i];
    Resistive_eta(v, grid->x[IDIR][i], grid->x[JDIR][j], grid->x[KDIR][k], J, eta);
    #endif

    // Poynting Flux
    Sx1[k][j][i] = (
      + d->Vc[VX3][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX3][k][j][i]
      - d->Vc[VX1][k][j][i] * pow(d->Vc[BX3][k][j][i], 2)
      - d->Vc[VX1][k][j][i] * pow(d->Vc[BX2][k][j][i], 2)
      + d->Vc[VX2][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX2][k][j][i]
      #if RESISTIVITY != NO
      + eta[IDIR] * (- Jx2[k][j][i] * d->Vc[BX3][k][j][i] + Jx3[k][j][i] * d->Vc[BX2][k][j][i])
      #endif
      );
    Sx2[k][j][i] = (
      + d->Vc[VX1][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX1][k][j][i]
      - d->Vc[VX2][k][j][i] * pow(d->Vc[BX1][k][j][i], 2)
      - d->Vc[VX2][k][j][i] * pow(d->Vc[BX3][k][j][i], 2)
      + d->Vc[VX3][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX3][k][j][i]
      #if RESISTIVITY != NO
      + eta[JDIR] * (- Jx3[k][j][i] * d->Vc[BX1][k][j][i] + Jx1[k][j][i] * d->Vc[BX3][k][j][i])
      #endif
      );
    Sx3[k][j][i] = (
      + d->Vc[VX2][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX2][k][j][i]
      - d->Vc[VX3][k][j][i] * pow(d->Vc[BX2][k][j][i], 2)
      - d->Vc[VX3][k][j][i] * pow(d->Vc[BX1][k][j][i], 2)
      + d->Vc[VX1][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX1][k][j][i]
      #if RESISTIVITY != NO
      + eta[KDIR] * (- Jx1[k][j][i] * d->Vc[BX2][k][j][i] + Jx2[k][j][i] * d->Vc[BX1][k][j][i])
      #endif
      );

    // Enthalpy flux
    E = (
      + 0.5 * d->Vc[RHO][k][j][i] * (pow(d->Vc[VX1][k][j][i], 2) + pow(d->Vc[VX2][k][j][i], 2) + pow(d->Vc[VX3][k][j][i], 2))
      + 2.5 * d->Vc[PRS][k][j][i]
      + d->Vc[RHO][k][j][i] * CONST_gsun * (2. * L / CONST_PI) * cos(0.5 * CONST_PI * grid->x[KDIR][k] / L)
      );
    Fx1[k][j][i] = - E * d->Vc[VX1][k][j][i];
    Fx2[k][j][i] = - E * d->Vc[VX2][k][j][i];
    Fx3[k][j][i] = - E * d->Vc[VX3][k][j][i];

  }
  */

}
/* ************************************************************* */
void ChangeOutputVar ()
/*
 *
 *
 *************************************************************** */
{
}
