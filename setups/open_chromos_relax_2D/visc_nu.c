/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
  double logRe;
  double tstop = RuntimeGet()->tstop;
  *nu1 = 1e8*1e7*UNIT_DENSITY/g_inputParam[Re]; // vscale*lscale*rhoscale/Re
  if (g_time > g_inputParam[t_Re_ramp]) {
    logRe = log(g_inputParam[Re])
      + log(g_inputParam[Re_end] / g_inputParam[Re])
      * (g_time - g_inputParam[t_Re_ramp]) / (tstop - g_inputParam[t_Re_ramp]);
    *nu1 = 1e8*1e7*UNIT_DENSITY/exp(logRe);
  }
  if (x2 < 60.) {
    *nu1 = 1e8*1e7*UNIT_DENSITY/1e8 * atan((60. - x2)/0.1) / CONST_PI * 2.;
  }
  *nu2 = 0.0;
}
