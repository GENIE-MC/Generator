/*
 * =====================================================================================
 *
 *       Filename:  onebody_currents_sf.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/09/2025 01:48:15 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIEa
 *
 *  \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
 *            For the full text of the license visit http://copyright.genie-mc.org
 *            or see $GENIE/LICENSE
 *
 * =====================================================================================
 */
#ifndef __ONEBODY_CURRENTS_SF_H__
#define __ONEBODY_CURRENTS_SF_H__
 
#include <complex>
#include <array>
#include <cmath>
 
namespace genie {
namespace onebody_currents_sf {
 
  void ComputeNucleonTensor(double xmn_in, double w_in, double wt,
                                double xk_x, double xk_y, double xk_z,
                                double q_x,  double q_y,  double q_z,
                                double ff1v, double ff2v,
                                double ffa,  double ffp,
                                std::complex<double> HadronTensor[4][4]);
 
  // current conservation
  void ComputeNucleonTensorCC(double xmn_in, double w_in, double wt,
                                  double xk_x, double xk_y, double xk_z,
                                  double q_x,  double q_y,  double q_z,
                                  double ff1v, double ff2v,
                                  double ffa,  double ffp,
                                  std::complex<double> HadronTensor[4][4]);
 
}  // namespace onebody_currents_sf
}  // namespace genie
 
#endif  // __ONEBODY_CURRENTS_SF_H__
