/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define CONST_gsun  (27.4e3 / (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH))

#define IDOM_LOOP_REV(i)  for ((i) = IEND; (i) >= IBEG; (i)--)
#define JDOM_LOOP_REV(j)  for ((j) = JEND; (j) >= JBEG; (j)--)
#define KDOM_LOOP_REV(k)  for ((k) = KEND; (k) >= KBEG; (k)--)

#define DOM_LOOP_REV(k,j,i) KDOM_LOOP_REV(k) JDOM_LOOP_REV(j) IDOM_LOOP_REV(i)

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  v[RHO] = 0.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
    v[PRS] = 0.0;
  #endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;

  #endif
}

/* ********************************************************************* */
double step_function(double x1, double x2, double r) {
/*!
 * Compute the step function that defines the loop section in the in the x1-x2
 * plane.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] r   ratio between the inside and the outside of the step
 * 
 * \return The value of the step function.
 *********************************************************************** */
  double step;
  step = 0.5*(1 - tanh( ( (sqrt( pow(x1,2) + pow(x2,2) )/g_inputParam[STEP_R]) - 1 ) * g_inputParam[STEP_b]));
  step = 1 + (r - 1) * step;
  return step;
}

void InitVariableFromDbl (Data *d, Grid *grid, char *data_dir, char *var, int nv) 
/*!
 * Initialize a given variable from dbl files
 *
 * \param [out] d         the PLUTO Data structure
 * \param [out] grid      pointer to array of Grid structures
 * \param [in] data_dir   path to the directory containing dbl and grid files
 * \param [in] var        name of the variable in the dbl files (eg. "prs")
 * \param [in] nv         index of the variable in d (eg. PRS)
 *
 *********************************************************************** */
{
  int k, j, i, dbl_fid;
  double *x1, *x2, *x3;
  char *data_fname, *grid_fname;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  char sep[] = "/";
  char data_fname_end[] = ".0000.dbl";
  char grid_fname_end[] = "grid.out";

  data_fname = malloc(strlen(data_dir) + strlen(sep) + strlen(var) + strlen(data_fname_end) + 1);
  strcpy(data_fname, data_dir);
  strcat(data_fname, sep);
  strcat(data_fname, var);
  strcat(data_fname, data_fname_end);
  grid_fname = malloc(strlen(data_dir) + strlen(sep) + strlen(grid_fname_end) + 1);
  strcpy(grid_fname, data_dir);
  strcat(grid_fname, sep);
  strcat(grid_fname, grid_fname_end);
  PlutoError(strlen(data_fname) > 63, "data_fname is too long");

  dbl_fid = InputDataOpen(data_fname, grid_fname, " ", 0);
  TOT_LOOP(k, j, i) {
    d->Vc[nv][k][j][i] = InputDataInterpolate(dbl_fid, x1[i], x2[j], x3[k]);
  }
  InputDataClose(dbl_fid);
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{

  char data_dir[] = "initial_state";
  InitVariableFromDbl(d, grid, data_dir, "prs", PRS);
  InitVariableFromDbl(d, grid, data_dir, "rho", RHO);
  // InitVariableFromDbl(d, grid, data_dir, "vx1", VX1);
  InitVariableFromDbl(d, grid, data_dir, "vx2", VX2);
  InitVariableFromDbl(d, grid, data_dir, "vx3", VX3);
  InitVariableFromDbl(d, grid, data_dir, "Bx1", BX1);
  InitVariableFromDbl(d, grid, data_dir, "Bx2", BX2);
  InitVariableFromDbl(d, grid, data_dir, "Bx3", BX3);
  InitVariableFromDbl(d, grid, data_dir, "psi_glm", PSI_GLM);

  // initial perturbation -----------------------------------------------------
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double v0 = g_inputParam[DRIVER_v0] / UNIT_VELOCITY;
  double L = g_inputParam[HALF_LOOP_L];
  double R = g_inputParam[STEP_R];
  double b = g_inputParam[STEP_b];
  double Rlim_sq = pow(R + 3.*R/b, 2);

  TOT_LOOP(k,j,i){
    if ( (pow(x1[i], 2) + pow(x2[j], 2)) <= Rlim_sq ) {
      d->Vc[VX1][k][j][i] = v0 * cos(CONST_PI*x3[k]/(2.*L));
    }
    else {
      d->Vc[VX1][k][j][i] = 0.;
    }
  }
  // --------------------------------------------------------------------------


}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double Temp, strat;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == X3_END){
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

        // Temperature and pressure
        #if CH_BOUNDARY_PRS_RHO == CH_STRATIFICATION
        Temp = d->Vc[PRS][k-1][j][i] / d->Vc[RHO][k-1][j][i];
        strat = CONST_gsun * sin(0.5 * CONST_PI * x3[k-1] / g_inputParam[HALF_LOOP_L]);
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k-1][j][i] + fabs(x3[k]-x3[k-1]) * strat * d->Vc[RHO][k-1][j][i];
        d->Vc[RHO][k][j][i] = d->Vc[PRS][k][j][i] / Temp;
        #endif
        #if CH_BOUNDARY_PRS_RHO == CH_OUTFLOW
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][JEND][i];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][JEND][i];
        #endif

        // Velocity
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[VX3][k][j][i] = 0.0;

        // Magnetic field
        #if CH_BOUNDARY_B == CH_EXTRAPOLATION
        d->Vc[BX1][k][j][i] =  (1.0/11.0)*(0.0 + 2.0*d->Vc[BX1][k-3][j][i] - 9.0*d->Vc[BX1][k-2][j][i] + 18.0*d->Vc[BX1][k-1][j][i]);
        d->Vc[BX2][k][j][i] =  (1.0/11.0)*(0.0 + 2.0*d->Vc[BX2][k-3][j][i] - 9.0*d->Vc[BX2][k-2][j][i] + 18.0*d->Vc[BX2][k-1][j][i]);
        d->Vc[BX3][k][j][i] =  (1.0/11.0)*(0.0 + 2.0*d->Vc[BX3][k-3][j][i] - 9.0*d->Vc[BX3][k-2][j][i] + 18.0*d->Vc[BX3][k-1][j][i]);
        #endif
        #if CH_BOUNDARY_B == CH_OUTFLOW
        d->Vc[BX1][k][j][i] = d->Vc[BX1][KEND][j][i];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][KEND][j][i];
        d->Vc[BX3][k][j][i] = d->Vc[BX2][KEND][j][i];
        #endif

      }
    }
  }

}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = CONST_gsun * sin( 0.5*CONST_PI*x3/g_inputParam[HALF_LOOP_L]);
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
