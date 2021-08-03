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

#define IN_SUBDOMAIN(x, i, beg, end) (x[i] >= beg && x[i] <= end)
#define ON_BEG_FACE(x, i, beg) (x[i-1] < beg && x[i] >= beg)
#define ON_END_FACE(x, i, end) (x[i] <= end && x[i+1] > end)


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
  InitVariableFromDbl(d, grid, data_dir, "vx1", VX1);
  InitVariableFromDbl(d, grid, data_dir, "vx2", VX2);
  InitVariableFromDbl(d, grid, data_dir, "vx3", VX3);
  InitVariableFromDbl(d, grid, data_dir, "Bx1", BX1);
  InitVariableFromDbl(d, grid, data_dir, "Bx2", BX2);
  InitVariableFromDbl(d, grid, data_dir, "Bx3", BX3);
  InitVariableFromDbl(d, grid, data_dir, "psi_glm", PSI_GLM);

}

/* ********************************************************************* */
double KineticEnergy(const double *v)
/*!
 *  Compute the kinetic energy.
 *
 * \param [in] v   a pointer to a vector of primitive variables
 * \return The kinetic energy \f$ 0.5 \rho \times (v_1^2+v_2^2+v_3^2) \f$.
 *
 *********************************************************************** */
{
  return 0.5*v[RHO] * (v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]);
}

double InternalEnergy(const double *v)
/*!
 *  Compute the internal energy.
 *
 * \param [in] v   a pointer to a vector of primitive variables
 * \return The internal energy \f$ p / (\gamma - 1) \f$.
 *
 *********************************************************************** */
{
  return v[PRS]/(g_gamma - 1.0);
}

double MagneticEnergy(const double *v)
/*!
 *  Compute the magnetic energy.
 *
 * \param [in] v   a pointer to a vector of primitive variables
 * \return The magnetic energy \f$ (B_1^2+B_2^2+B_3^2)/2 \f$.
 *
 *********************************************************************** */
{
  return 0.5 * (v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
}

double GravitationalEnergy(const double *v, double x1, double x2, double x3)
/*!
 *  Compute the gravitational energy.
 *
 * \param [in] v   a pointer to a vector of primitive variables
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \return The gravitational energy \f$ \rho \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  #if BODY_FORCE == NO
  return 0.;
  #elif BODY_FORCE == POTENTIAL
  return v[RHO] * BodyForcePotential(x1, x2, x3);
  #else
  printf(" ! Can't compute gravitational energy with non-potential body force\n");
  QUIT_PLUTO(1);
  #endif
}

void ElectricCurrent(const Data *d, Grid *grid, double *J, int i, int j, int k)
/*!
 *  Compute the electric current.
 *
 * \param [in] d  the PLUTO Data structure
 * \param [in] grid  pointer to an array of Grid structures.
 * \param [out] J the electric current
 * \param [in] i  index in the 1st coordinate direction
 * \param [in] j  index in the 2nd coordinate direction
 * \param [in] k  index in the 3rd coordinate direction
 *
 *********************************************************************** */
{
  J[IDIR] = 0.5*(
    + (d->Vc[BX3][k][j+1][i] - d->Vc[BX3][k][j-1][i]) / grid->dx[JDIR][j]
    - (d->Vc[BX2][k+1][j][i] - d->Vc[BX2][k-1][j][i]) / grid->dx[KDIR][k]);
  J[JDIR] = 0.5*(
    - (d->Vc[BX3][k][j][i+1] - d->Vc[BX3][k][j][i-1]) / grid->dx[IDIR][i]
    + (d->Vc[BX1][k+1][j][i] - d->Vc[BX1][k-1][j][i]) / grid->dx[KDIR][k]);
  J[KDIR] = 0.5*(
    + (d->Vc[BX2][k][j][i+1] - d->Vc[BX2][k][j][i-1]) / grid->dx[IDIR][i]
    - (d->Vc[BX1][k][j+1][i] - d->Vc[BX1][k][j-1][i]) / grid->dx[JDIR][j]);
}

void PoyntingFlux(const Data *d, Grid *grid, double *S, int i, int j, int k)
/*!
 *  Compute the Poynting flux.
 *
 * \param [in] d  the PLUTO Data structure
 * \param [in] grid  pointer to an array of Grid structures.
 * \param [out] S the Poynting flux
 * \param [in] i  index in the 1st coordinate direction
 * \param [in] j  index in the 2nd coordinate direction
 * \param [in] k  index in the 3rd coordinate direction
 *
 *********************************************************************** */
{
  // Resistivity
  #if RESISTIVITY != NO
  int nv;
  double *v;
  double J[3], eta[3];
  v = ARRAY_1D(NVAR, double);
  ElectricCurrent(d, grid, J, i, j, k);
  NVAR_LOOP(nv) { v[nv] = d->Vc[nv][k][j][i]; }
  Resistive_eta(v, grid->x[IDIR][i], grid->x[JDIR][j], grid->x[KDIR][k], J, eta);
  #endif

  S[IDIR] += (
    + d->Vc[VX3][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX3][k][j][i]
    - d->Vc[VX1][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX3][k][j][i]
    - d->Vc[VX1][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX2][k][j][i]
    + d->Vc[VX2][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX2][k][j][i]
    #if RESISTIVITY != NO
    + eta[IDIR] * (- J[JDIR] * d->Vc[BX3][k][j][i] + J[KDIR] * d->Vc[BX2][k][j][i])
    #endif
    );

  S[JDIR] = (
    + d->Vc[VX1][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX1][k][j][i]
    - d->Vc[VX2][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX1][k][j][i]
    - d->Vc[VX2][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX3][k][j][i]
    + d->Vc[VX3][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX3][k][j][i]
    #if RESISTIVITY != NO
    + eta[JDIR] * (- J[KDIR] * d->Vc[BX1][k][j][i] + J[IDIR] * d->Vc[BX3][k][j][i])
    #endif
    );

  S[KDIR] = (
    + d->Vc[VX2][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX2][k][j][i]
    - d->Vc[VX3][k][j][i] * d->Vc[BX2][k][j][i] * d->Vc[BX2][k][j][i]
    - d->Vc[VX3][k][j][i] * d->Vc[BX1][k][j][i] * d->Vc[BX1][k][j][i]
    + d->Vc[VX1][k][j][i] * d->Vc[BX3][k][j][i] * d->Vc[BX1][k][j][i]
    #if RESISTIVITY != NO
    + eta[KDIR] * (- J[IDIR] * d->Vc[BX2][k][j][i] + J[JDIR] * d->Vc[BX1][k][j][i])
    #endif
    );

}

void EnthalpyFlux(const double *v, double *F, double x1, double x2, double x3)
/*!
 *  Compute the enthalpy flux.
 *
 * \param [in] v   a pointer to a vector of primitive variables
 * \param [out] F the enthalpy flux
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  double E;
  E = (KineticEnergy(v)
       + g_gamma * InternalEnergy(v)
       + GravitationalEnergy(v, x1, x2, x3));
  F[IDIR] = - E * v[VX1];
  F[JDIR] = - E * v[VX2];
  F[KDIR] = - E * v[VX3];
}

int FluxIndex(int dir, int side)
{
  return dir*2 + side;
}

void EnergyAnalysis (const Data *d, Grid *grid, char* output_file, double x1beg, double x1end, double x2beg, double x2end, double x3beg, double x3end)
/*!
 *  Perform runtime energy analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 * \param [in] x1beg, x1end, x2beg, x2end, x3beg, x3end  boundaries of the
 * subdomain (use NAN to get compute over the full domain).
 *
 *********************************************************************** */
{

  int i, j, k, nv, dir, side;
  double *v;
  double Lx, Ly, Lz, dA, dV, vol;
  double scrhv[3];
  double EK, EI, EB, EG;
  double S[6], F[6];
  double send[16], recv[16];
  double *x, *y, *z;
  double *dx, *dy, *dz;

  v = ARRAY_1D(NVAR, double);

  /* ---- Set pointer shortcuts ---- */
  x = grid->x[IDIR];
  y = grid->x[JDIR];
  z = grid->x[KDIR];
  dx = grid->dx[IDIR];
  dy = grid->dx[JDIR];
  dz = grid->dx[KDIR];

  /* ---- Sub domain bounds ---- */

  if (isnan(x1beg)) x1beg = g_domBeg[IDIR];
  if (isnan(x1end)) x1end = g_domEnd[IDIR];
  if (isnan(x2beg)) x2beg = g_domBeg[JDIR];
  if (isnan(x2end)) x2end = g_domEnd[JDIR];
  if (isnan(x3beg)) x3beg = g_domBeg[KDIR];
  if (isnan(x3end)) x3end = g_domEnd[KDIR];

  /* ---- Main loops ---- */
  EK = 0.;
  EI = 0.;
  EB = 0.;
  EG = 0.;
  for (dir = 0; dir < 3; dir++) {
    for (side = 0; side < 2; side++) {
      S[FluxIndex(dir, side)] = 0.;
      F[FluxIndex(dir, side)] = 0.;
    }
  }

  // volume
  DOM_LOOP(k,j,i) {
    if (IN_SUBDOMAIN(x, i, x1beg, x1end) &&
        IN_SUBDOMAIN(y, j, x2beg, x2end) &&
        IN_SUBDOMAIN(z, k, x3beg, x3end)) {
      NVAR_LOOP(nv) { v[nv] = d->Vc[nv][k][j][i]; }

      dV = dx[i]*dy[j]*dz[k]; /* Cell volume (Cartesian coordinates) */

      EK += dV * KineticEnergy(v);
      EI += dV * InternalEnergy(v);
      EB += dV * MagneticEnergy(v);
      EG += dV * GravitationalEnergy(v, x[i], y[j], z[k]);
    }
  }

  // X3
  DOM_LOOP(k,j,i) {
    if (IN_SUBDOMAIN(x, i, x1beg, x1end) &&
        IN_SUBDOMAIN(y, j, x2beg, x2end)) {

      if      ON_BEG_FACE(z, k, x3beg) side = 0;
      else if ON_END_FACE(z, k, x3end) side = 1;
      else side = -1;

      if (side >= 0) {
        dA = dx[i]*dy[j]; /* Cell surface */
        NVAR_LOOP(nv) { v[nv] = d->Vc[nv][k][j][i]; }
        PoyntingFlux(d, grid, scrhv, i, j, k);
        S[FluxIndex(KDIR, side)] += dA * scrhv[KDIR];
        EnthalpyFlux(v, scrhv, x[i], y[j], z[k]);
        F[FluxIndex(KDIR, side)] += dA * scrhv[KDIR];
      }

    }
  }

  // X2
  DOM_LOOP(k,j,i) {
    if (IN_SUBDOMAIN(x, i, x1beg, x1end) &&
        IN_SUBDOMAIN(z, k, x3beg, x3end)) {

      if      ON_BEG_FACE(y, j, x2beg) side = 0;
      else if ON_END_FACE(y, j, x2end) side = 1;
      else side = -1;

      if (side >= 0) {
        dA = dx[i]*dz[k]; /* Cell surface */
        NVAR_LOOP(nv) { v[nv] = d->Vc[nv][k][j][i]; }
        PoyntingFlux(d, grid, scrhv, i, j, k);
        S[FluxIndex(JDIR, side)] += dA * scrhv[JDIR];
        EnthalpyFlux(v, scrhv, x[i], y[j], z[k]);
        F[FluxIndex(JDIR, side)] += dA * scrhv[JDIR];
      }

    }
  }

  // X1
  DOM_LOOP(k,j,i) {
    if (IN_SUBDOMAIN(y, j, x2beg, x2end) &&
        IN_SUBDOMAIN(z, k, x3beg, x3end)) {

      if      ON_BEG_FACE(x, i, x1beg) side = 0;
      else if ON_END_FACE(x, i, x1end) side = 1;
      else side = -1;

      if (side >= 0) {
        dA = dy[j]*dz[k]; /* Cell surface */
        NVAR_LOOP(nv) { v[nv] = d->Vc[nv][k][j][i]; }
        PoyntingFlux(d, grid, scrhv, i, j, k);
        S[FluxIndex(IDIR, side)] += dA * scrhv[IDIR];
        EnthalpyFlux(v, scrhv, x[i], y[j], z[k]);
        F[FluxIndex(IDIR, side)] += dA * scrhv[IDIR];
      }

    }
  }

  /* Compute total domain volume and surfaces */
  Lx = x1end - x1beg,
  Ly = x2end - x2beg;
  Lz = x3end - x3beg;
  vol = Lx * Ly * Lz;

  // Compute energy average
  EK /= vol;
  EI /= vol;
  EB /= vol;
  EG /= vol;
  for (side = 0; side < 2; side++) {
    S[FluxIndex(IDIR, side)] /= vol;
    S[FluxIndex(JDIR, side)] /= vol;
    S[FluxIndex(KDIR, side)] /= vol;
    F[FluxIndex(IDIR, side)] /= vol;
    F[FluxIndex(JDIR, side)] /= vol;
    F[FluxIndex(KDIR, side)] /= vol;
  }

  // Compute flux average

  /* ---- Parallel data reduction ---- */
  #ifdef PARALLEL
    send[0] = EK;
    send[1] = EI;
    send[2] = EB;
    send[3] = EG;
    send[4+0] = S[0];
    send[4+1] = S[1];
    send[4+2] = S[2];
    send[4+3] = S[3];
    send[4+4] = S[4];
    send[4+5] = S[5];
    send[10+0] = F[0];
    send[10+1] = F[1];
    send[10+2] = F[2];
    send[10+3] = F[3];
    send[10+4] = F[4];
    send[10+5] = F[5];

    MPI_Allreduce(&send, &recv, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    EK = recv[0];
    EI = recv[1];
    EB = recv[2];
    EG = recv[3];
    S[0] = recv[4+0];
    S[1] = recv[4+1];
    S[2] = recv[4+2];
    S[3] = recv[4+3];
    S[4] = recv[4+4];
    S[5] = recv[4+5];
    F[0] = recv[10+0];
    F[1] = recv[10+1];
    F[2] = recv[10+2];
    F[3] = recv[10+3];
    F[4] = recv[10+4];
    F[5] = recv[10+5];

    MPI_Barrier (MPI_COMM_WORLD);
  #endif

  /* ---- Write ascii file to disk ---- */
  if (prank == 0) {
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;

    sprintf(fname, "%s/%s",RuntimeGet()->output_dir,output_file);
    if (g_stepNumber == 0) { /* Open for writing only when we’re starting */
      fp = fopen(fname,"w"); /* from beginning */
      fprintf (fp,"# %10s %12s %22s %22s %22s %22s","t","dt","EB","EI","EK","EG");
      fprintf(fp," %22s %22s %22s %22s %22s %22s","SxBEG","SxEND","SyBEG","SyEND","SzBEG","SzEND");
      fprintf(fp," %22s %22s %22s %22s %22s %22s","FxBEG","FxEND","FyBEG","FyEND","FzBEG","FzEND");
      fprintf(fp,"\n");
      /* Write log */
      print("> Energy analysis:\n");
      print("  Output file: %s\n", fname);
      print("  Domain:\n");
      print("  X1: [ %8.4f, %8.4f]\n", x1beg, x1end);
      print("  X2: [ %8.4f, %8.4f]\n", x2beg, x2end);
      print("  X3: [ %8.4f, %8.4f]\n", x3beg, x3end);
      print("\n");
    } else {
      /* Append if this is not step 0 */
      if (tpos < 0.0) { /* Obtain time coordinate of to last written row */
        char sline[512];
        fp = fopen(fname,"r");
        if (fp != NULL) {
          while (fgets(sline, 512, fp)) {}
          sscanf(sline, "%lf\n",&tpos); /* tpos = time of the last written row */
          fclose(fp);
        }
      }
      fp = fopen(fname,"a");
    }
    if (g_time > tpos) { /* Write if current time if > tpos */
      fprintf (fp,"%.6e %.6e %+.15e %+.15e %+.15e %+.15e",g_time,g_dt,EB,EI,EK,EG);
      fprintf(fp," %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e",S[0],S[1],S[2],S[3],S[4],S[5]);
      fprintf(fp," %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e",F[0],F[1],F[2],F[3],F[4],F[5]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

}

void CutZAnalysis(const Data *d, Grid *grid, char* output_file, double x1cut, double x2cut)
/*!
 *  Save cut along a lign
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 * \param [in] x1, x2  coordinates of the cut
 * subdomain (use NAN to get compute over the full domain).
 *
 *********************************************************************** */
{
  DEBUG_FUNC_BEG("CutZAnalysis");

  int i, j, k, nv;

  // -- Pointer shortcuts
  double *x, *y;
  x = grid->x[IDIR];
  y = grid->x[JDIR];

  // -- Determine cut line goes through this process' subdomain
  int contains_cut, contains_x1cut, contains_x2cut;
  contains_x1cut = ((grid->xbeg[IDIR] < x1cut) && (x1cut <= grid->xend[IDIR]));
  contains_x2cut = ((grid->xbeg[JDIR] < x2cut) && (x2cut <= grid->xend[JDIR]));
  contains_cut = contains_x1cut & contains_x2cut;

  // -- Determine cut indices
  int i0 = -1;
  int j0 = -1;
  double delta;
  if (contains_cut) {
    delta = 1000.;
    IDOM_LOOP(i) {
      if (fabs(x[i] - x1cut) < delta) {
        delta = fabs(x[i] - x1cut);
        i0 = i;
      }
    }
    delta = 1000.;
    JDOM_LOOP(j) {
      if (fabs(y[j] - x2cut) < delta) {
        delta = fabs(y[j] - x2cut);
        j0 = j;
      }
    }
  }

  #if DEBUG == TRUE
  print("    [CutZ slice (%+.6g %+.6g %+.6g)>%d (%.6g %.6g %.6g)>%d (%d %d) (%.6g %.6g)]\n",
    grid->xbeg[IDIR], x1cut, grid->xend[IDIR], contains_x1cut,
    grid->xbeg[JDIR], x2cut, grid->xend[JDIR], contains_x2cut,
    i0, j0,
    x[i0], y[j0]);
  #endif

  // -- Cut array size (local)
  int n_VcZ = grid->np_int[KDIR];

  // -- Init receiving cut array (root process)
  double *VcZ_glob;
  int n_VcZ_glob = grid->np_int_glob[KDIR];
  if (prank == 0) {
    VcZ_glob = ARRAY_1D(NVAR * n_VcZ_glob, double);
  }
  else {
    VcZ_glob = NULL;
  }

  // -- Gather data on root process

  // ---- Data counts (for MPI_Gatherv)
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  // ------ Sent data
  // Only data if cut line goes through the process' subdomain
  double *sendbuf;
  int sendcount;
  if (contains_cut) {
    sendbuf = ARRAY_1D(n_VcZ, double);
    sendcount = n_VcZ;
  }
  else {
    sendbuf = NULL;
    sendcount = 0;
  }
  // ------ Received data
  // Receive data with root process
  double *recvbuf;
  int *recvcounts, *displs;
  if (prank == 0) {
    recvbuf = ARRAY_1D(n_VcZ_glob, double);
    recvcounts = ARRAY_1D(nproc, int);
    displs = ARRAY_1D(nproc, int);
  }
  else {
    recvbuf = NULL;
    recvcounts = NULL;
    displs = NULL;
  }
  // gather sent data size on root process (needed by MPI_Gatherv)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&sendcount, 1, MPI_INT, (void *)recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  // determine displacements (for gatherv)
  if (prank == 0) {
    int offset = 0;
    int n;
    for (n = 0; n < nproc; n++) {
      displs[n] = offset;
      offset += recvcounts[n];
    }
  }

  // ---- Exchange data
  NVAR_LOOP(nv) {
    if (contains_cut) {
      KDOM_LOOP(k) {
        sendbuf[k-KBEG] = d->Vc[nv][k][j0][i0];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv((void *)sendbuf, sendcount, MPI_DOUBLE, (void *)recvbuf, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (prank == 0) {
      for (k = 0; k < n_VcZ_glob; k++) {
        VcZ_glob[k + n_VcZ_glob * nv] = recvbuf[k];
      }
    }
  }


  // -- Write data to disk
  if (prank == 0) {
    char filename[512];
    FILE *fp_list, *fp_bin;

    // ---- Open list file
    double tpos = -1.;
    int nfile = -1;
    sprintf(filename, "%s/%s.list.out", RuntimeGet()->output_dir, output_file);
    #if DEBUG == TRUE
      print("    [CutZ output %s ", filename);
    #endif
    if (g_stepNumber == 0) {
      // Open file for writing when starting sim from 0
      #if DEBUG == TRUE
        print("A ");
      #endif
      fp_list = fopen(filename, "w");
    } else {
      // Append to file if not starting sim from 0.
      // In this case, time coordinate of to last written row when starting the
      // simulation.
      #if DEBUG == TRUE
        print("B ");
      #endif
      if (tpos < 0) {
        char sline[512];
        fp_list = fopen(filename, "r");
        if (fp_list != NULL) {
          while (fgets(sline, 512, fp_list)) {}
          sscanf(sline, "%d %lf\n", &nfile, &tpos); // tpos = time of the last written row
          fclose(fp_list);
        }
        #if DEBUG == TRUE
          print("nfile=%d tpos=%g ", nfile, tpos);
        #endif
      }
      fp_list = fopen(filename, "a");
    }
    #if DEBUG == TRUE
      print("]\n");
    #endif

    // ---- Write data
    if (g_time > tpos) {
      // Write if current time if > tpos
      nfile += 1;

      // ----- Write to list file
      fprintf(fp_list, "%d %16.10e %12.6e %ld\n", nfile, g_time, g_dt, g_stepNumber);

      // ----- Write binary data
      // ------ Determine output filename
      sprintf(filename, "%s/%s.%04d.dat", RuntimeGet()->output_dir, output_file, nfile);
      // ------ Write data
      fp_bin = fopen(filename, "wb");
      fwrite(VcZ_glob, sizeof(double), NVAR * n_VcZ_glob, fp_bin);
      fclose(fp_bin);

      // ----- log message
      print("> Writing CutZ file %s to disk (x=%g y=%g t=%g)]\n",
        filename, x[i0], y[j0], g_time);

    }

    // ---- Close list file
    fclose(fp_list);

  }

  DEBUG_FUNC_END("CutZAnalysis");
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
  EnergyAnalysis(d, grid, "energy_all.dat", NAN, NAN, NAN, NAN, NAN, NAN);
  EnergyAnalysis(d, grid, "energy_cor.dat", NAN, NAN, NAN, NAN, NAN, 90.);
  CutZAnalysis(d, grid, "cut_center_pp", +.02, +.02);
  CutZAnalysis(d, grid, "cut_center_pm", +.02, -.02);
  CutZAnalysis(d, grid, "cut_center_mp", -.02, +.02);
  CutZAnalysis(d, grid, "cut_center_mm", -.02, -.02);
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
  // retrieve input parameters
  const double av_layer_min = g_inputParam[VRW_AV_LAYER_MIN];
  const double x_max_layer = g_inputParam[VRW_X_MAX_LAYER];

  // bv boundary layer (x3 dependence)
  double bv_layer;
  if (x_loop < x_max_layer) {
    bv_layer = av_layer_min + (1. - av_layer_min ) / x_max_layer * x_loop;
  }
  else {
    bv_layer = 1.;
  }

  return bv_layer;
}

double* LoadTabVelocity(char *fname, int *n)
/*!
 *  Load tabulated driver velocity
 *
 * \param [in] fname  the name of the input data file
 * \param [out] n  number of velocity points returned
 *
 * \return Array containing the velocity at each timestep, normalized to its
 *         variance
 *
 *********************************************************************** */
{

  FILE *fp;
  char sline[512];

  // Get number of lines
  (*n) = 0;
  fp = fopen(fname, "r");
  if (fp != NULL) {
    while (fgets(sline, 512, fp)) {
      (*n) += 1;
    }
    fclose(fp);
  }
  else {
    print("! Can't read file %s\n", fname);
    QUIT_PLUTO(1);
  }

  // Allocate memory
  double *v = ARRAY_1D((*n), double);

  // Load data
  double v_in;
  fp = fopen(fname, "r");
  if (fp != NULL) {
    for (int i = 0; i < (*n); i++) {
      // I assume that the number of lines (*n) won't change since counting
      // them, because the velocity file is only an input. If it wasn't the
      // case, this could trigger some segfault.
      fgets(sline, 512, fp);
      sscanf(sline, "%lf\n", &v_in);
      v[i] = v_in;
    }
    fclose(fp);
  }
  else {
    print("! Can't read file %s\n", fname);
    QUIT_PLUTO(1);
  }

  return v;
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
  int   i, j, k;
  double  *x1, *x2, *x3;
  double vnew, x1new, x2new, profile, profile1, profile2, Temp, strat;

  double radius = g_inputParam[STEP_R];
  double v0_driver = g_inputParam[DRIVER_v0] / UNIT_VELOCITY;

  // Setup broadband driver
  // load data
  static double *tab_vnew = NULL;
  static int n_vnew;
  if (tab_vnew == NULL) {
    tab_vnew = LoadTabVelocity("driver_v/v.txt", &n_vnew);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // velocity amplitude at current timestep
  vnew = v0_driver * tab_vnew[g_stepNumber];


  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  /* Lower boundary */
  if (side == X3_END){
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

        // Temperature and pressure
        #if CH_BOUNDARY_PRS_RHO == CH_STRATIFICATION
        Temp = d->Vc[PRS][k-1][j][i] / d->Vc[RHO][k-1][j][i];
        strat = CONST_gsun;
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k-1][j][i] + fabs(x3[k]-x3[k-1]) * strat * d->Vc[RHO][k-1][j][i];
        d->Vc[RHO][k][j][i] = d->Vc[PRS][k][j][i] / Temp;
        #endif
        #if CH_BOUNDARY_PRS_RHO == CH_OUTFLOW
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][JEND][i];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][JEND][i];
        #endif

        // Velocity (driver)
        x1new = x1[i];  // + coordinate change
        x2new = x2[j];
        profile = 0.5 * (1 - tanh(((sqrt(pow(x1new, 2) + pow(x2new, 2)) / radius) - 1 ) * g_inputParam[STEP_b]));
        profile1 = pow(radius, 2) * (pow(x1new, 2) - pow(x2new, 2)) / pow(pow(x1new, 2) + pow(x2new, 2), 2);
        profile2 = pow(radius, 2) * 2 * x1new * x2new / pow(pow(x1new, 2) + pow(x2new, 2), 2);

        d->Vc[VX1][k][j][i] = vnew * (profile1 + (1.0 - profile1) * profile);
        d->Vc[VX2][k][j][i] = vnew * (profile2 + (0.0 - profile2) * profile);
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

  /* Velocity rewrite in the entire domain */
  if (side == 0){
    TOT_LOOP(k, j, i) {
      double av = VelocityRewriteCoeff(x3[k]);
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
  g[JDIR] = 0.0;
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
  double L = g_inputParam[HALF_LOOP_L];
  return CONST_gsun * (L - x3);
}
#endif

// vim: fdm=syntax
