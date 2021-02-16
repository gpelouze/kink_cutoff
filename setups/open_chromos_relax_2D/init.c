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

#define input_Bz0  (g_inputParam[Bz0]/sqrt(4*CONST_PI*UNIT_DENSITY)/UNIT_VELOCITY)
#define RHO_CH  (g_inputParam[NE_CH] * CONST_amu / UNIT_DENSITY)

#define pulseAmp  6e5/UNIT_VELOCITY
#define period  24.71

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
    v[BX2] = input_Bz0;
    v[BX3] = 0.0;

  #endif
}

/* ********************************************************************* */
double step_function(double x1, double r) {
/*!
 * Compute the step function that defines the loop section along the x1
 * direction.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] r   ratio between the inside and the outside of the step
 * 
 * \return The value of the step function.
 *********************************************************************** */
  double step;
  step = 0.5*(1 - tanh( ( (fabs(x1)/g_inputParam[STEP_R]) - 1 ) * g_inputParam[STEP_b]));
  step = 1 + (r - 1) * step;
  return step;
}

/* ********************************************************************* */
void HsEquilBodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the body force vector to be used by InitDomainSolveHsEquil from
 * vector and potential body force definitions.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  #if BODY_FORCE == NO
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
  #elif BODY_FORCE == VECTOR
  BodyForceVector(v, g, x1, x2, x3);
  #else
  printf(" ! Initial hs equil is not implemented for non-vector body force\n");
  QUIT_PLUTO(1);
  #endif
}

/* ********************************************************************* */
void InitDomainSolveHsEquil (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i, j, k, nv, dbl_fid;
  double *v, *x1, *x2, *x3;
  double g[3];
  double ***T, ***logPrs;
  double *T_1D;
  double step_T, step_prs, prs_ch, Tv, v1, prs, ds, mu;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  v = ARRAY_1D(NVAR, double);
  T = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  T_1D = ARRAY_1D(NX2_TOT, double);
  logPrs = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  // // Load and interpolate the temperature from dbl file
  dbl_fid = InputDataOpen("init_t_profile/T.0000.dbl", "init_t_profile/grid.out", "little", 0);
  // In the input 1D dbl file, the loop is along IDIR, but in this 2D setup,
  // the loop is along JDIR. Hence, InputDataInterpolate is called with the
  // x2 of this setup as the x1 of the input dbl file.
  JTOT_LOOP(j) {
    T_1D[j] = InputDataInterpolate(dbl_fid, x2[j], 0., 0.);
  }
  InputDataClose(dbl_fid);

  // Initialize T over full domain such that it is constant in the chromosphere
  // and stepped in the corona
  TOT_LOOP(k,j,i){
    step_T = step_function(x1[i], g_inputParam[STEP_rT]);
    T[k][j][i] = g_inputParam[T_CH] / KELVIN
      + (g_inputParam[T0_APEX] * step_T - g_inputParam[T_CH]) / KELVIN
      * T_1D[j];
  }

  // Solve hydrostatic equilibrium for pressure assuming constant B along x2
  // boundary condition at JEND
  KDOM_LOOP(k) IDOM_LOOP(i) {
    NVAR_LOOP(nv) {
      v[nv] = d->Vc[nv][k][JEND][i];
      HsEquilBodyForceVector(v, g, x1[i], x2[JEND], x3[k]);
      }
    step_prs = step_function(x1[i], g_inputParam[STEP_rprs]);
    mu = MeanMolecularWeight(v);
    ds = x2[JEND-1] - x2[JEND];
    Tv = (T[k][JEND][i] + T[k][JEND-1][i]) / 2.;
    v1 = ds * (- mu*g[JDIR]) / Tv;
    prs_ch = RHO_CH * T[k][JEND][i] / mu * step_prs;
    logPrs[k][JEND][i] = log(2. * prs_ch / (1. + exp(-v1)));
  }
  // solve along JDIR from JEND to JBEG
  DOM_LOOP_REV(k, j, i) {
    NVAR_LOOP(nv) {
      v[nv] = d->Vc[nv][k][j][i];
      HsEquilBodyForceVector(v, g, x1[i], x2[j], x3[k]);
      }
    mu = MeanMolecularWeight(v);
    ds = x2[j-1] - x2[j];
    Tv = (T[k][j][i] + T[k][j-1][i]) / 2.;
    logPrs[k][j-1][i] = logPrs[k][j][i] - ds * (- mu*g[JDIR]) / Tv;
  }

  // set rho and prs
  DOM_LOOP(k, j, i) {
    NVAR_LOOP(nv) {
      v[nv] = d->Vc[nv][k][j][i];
      }
    mu = MeanMolecularWeight(v);
    prs = exp(logPrs[k][j][i]);
    d->Vc[RHO][k][j][i] = prs / T[k][j][i] * mu;
    #if HAVE_ENERGY
      d->Vc[PRS][k][j][i] = prs;
    #endif
  }

}

/* ********************************************************************* */
void InitDomainLoadHsEquil (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i, j, k, dbl_fid;
  double *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  // Load and interpolate pressure and temperature from dbl files in
  // init_hs_equil/

  dbl_fid = InputDataOpen("init_hs_equil/rho.0000.dbl", "init_hs_equil/grid.out", " ", 0);
  TOT_LOOP(k, j, i) {
    d->Vc[RHO][k][j][i] = InputDataInterpolate(dbl_fid, x1[i], x2[j], x3[k]);
  }
  InputDataClose(dbl_fid);

  #if HAVE_ENERGY
  dbl_fid = InputDataOpen("init_hs_equil/prs.0000.dbl", "init_hs_equil/grid.out", " ", 0);
  TOT_LOOP(k, j, i) {
    d->Vc[PRS][k][j][i] = InputDataInterpolate(dbl_fid, x1[i], x2[j], x3[k]);
  }
  InputDataClose(dbl_fid);
  #endif

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
#if INIT_DOMAIN_FUNCTION == INIT_DOMAIN_SOLVE_HS_EQUIL
InitDomainSolveHsEquil(d, grid);
#endif
#if INIT_DOMAIN_FUNCTION == INIT_DOMAIN_LOAD_HS_EQUIL
InitDomainLoadHsEquil(d, grid);
#endif
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
double VelocityRewriteCoeff(double x_loop)
/*!
 *  Determine the velocity rewrite coefficient alpha_v.
 *
 * \param [in] x_loop  coordinate along the loop (x2 in 2D setups,
 *                     and x3 in 3D setups)
 *
 * \return the velocity rewrite coefficient \f$ \alpha_v(x_\mathrm{loop} \f$.
 *
 *********************************************************************** */
{
  // (beg) user-defined parameters
  // TODO: move to pluto.ini
  const double vrv_dt1 = 300.;  // hydrostatic equilibrium
  const double vrv_dt2 = 100.;
  const double vrv_dt3 = 300.;  // wave damping
  const double vrv_dt4 = 100.;
  const double vrv_dt5 = 300.;  // relaxation of damping
  const double vrv_dt6 = 100.;
  const double vrv_av_full_domain = .999999;
  const double vrv_av_min_boundary_layer = 0.;
  const double vrv_x_max_boundary_layer = 50.;
  // (end) user-defined parameters

  const double t1 = vrv_dt1;
  const double t2 = t1 + vrv_dt2;
  const double t3 = t2 + vrv_dt3;
  const double t4 = t3 + vrv_dt4;
  const double t5 = t4 + vrv_dt5;
  const double t6 = t5 + vrv_dt6;

  // av full domain (time dependence only)
  double av_relax;
  if (g_time < t1) {  // dt1
    av_relax = 1.;
  }
  else if (g_time < t2) {  // dt2
    av_relax = 1. - (g_time - t1) * (1. - vrv_av_full_domain) / vrv_dt2;
  }
  else if (g_time < t3) {  // dt3
    av_relax = vrv_av_full_domain;
  }
  else if (g_time < t4) {  // dt4
    av_relax = vrv_av_full_domain + (g_time - t3) * (1. - vrv_av_full_domain) / vrv_dt4;
  }
  else if (g_time < t5) {  // dt5
    av_relax = 1.;
  }
  else if (g_time < t6) {  // dt6
    av_relax = 1. - (g_time - t5) / vrv_dt6;
  }
  else {  // after dt6
    av_relax = 0.;
  }

  // av boundary layer (time dependence)
  double av_layer;
  if (g_time < t5) {  // dt1 to dt5
    av_layer = 0.;
  }
  else if (g_time < t6) {  // dt6
    av_layer = (g_time - t5) / vrv_dt6;
  }
  else {  // after dt6
    av_layer = 1.;
  }

  // av boundary layer (x2 dependence)
  double bv_layer;
  if (x_loop < vrv_x_max_boundary_layer) {
    bv_layer = vrv_av_min_boundary_layer + (1. - vrv_av_min_boundary_layer) / vrv_x_max_boundary_layer * x_loop;
  }
  else {
    bv_layer = 1.;
  }

  return av_relax + av_layer * bv_layer;
}

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

  double radius = g_inputParam[STEP_R];

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  /* Lower boundary */
  if (side == X2_END){
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

        // corrections for temperature-pressure
        #if CH_BOUNDARY_PRS_RHO == CH_STRATIFICATION
        Temp = d->Vc[PRS][k][j-1][i]/d->Vc[RHO][k][j-1][i];
        strat = CONST_gsun;
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j-1][i] + fabs(x2[j]-x2[j-1])*strat*d->Vc[RHO][k][j-1][i];
        d->Vc[RHO][k][j][i] =  d->Vc[PRS][k][j][i]/Temp;
        #endif
        #if CH_BOUNDARY_PRS_RHO == CH_OUTFLOW
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][JEND][i];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][JEND][i];
        #endif

        // velocity
        d->Vc[VX1][k][j][i] = + d->Vc[VX1][k][2*JEND-j+1][i];
        d->Vc[VX2][k][j][i] = - d->Vc[VX2][k][2*JEND-j+1][i];

        // for magnetic field extrapolation
        #if CH_BOUNDARY_B == CH_EXTRAPOLATION
        d->Vc[BX1][k][j][i] =  (1.0/11.0)*(0.0 + 2.0*d->Vc[BX1][k][j-3][i] - 9.0*d->Vc[BX1][k][j-2][i] + 18.0*d->Vc[BX1][k][j-1][i]);
        d->Vc[BX2][k][j][i] =  (1.0/11.0)*(0.0 + 2.0*d->Vc[BX2][k][j-3][i] - 9.0*d->Vc[BX2][k][j-2][i] + 18.0*d->Vc[BX2][k][j-1][i]);
        #endif
        #if CH_BOUNDARY_B == CH_OUTFLOW
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][JEND][i];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][JEND][i];
        #endif

      }
    }
  }

  /* Velocity rewrite in the entire domain */
  if (side == 0){
    TOT_LOOP(k, j, i) {
      double av = VelocityRewriteCoeff(x2[j]);
      EXPAND( d->Vc[VX1][k][j][i] *= av; ,
              d->Vc[VX2][k][j][i] *= av; ,
              d->Vc[VX3][k][j][i] *= av; )
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
  g[JDIR] = CONST_gsun;
  g[KDIR] = 0.0;
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
