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
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */

#ifndef __ONEBODY_CURRENTS_SF_H__
#define __ONEBODY_CURRENTS_SF_H__


#include <Eigen/Dense>
#include <complex>
#include <array>
#include <cmath>

namespace genie{
  namespace onebody_currents_sf {

    using cdouble = std::complex<double>;
    using Mat2c   = Eigen::Matrix2cd;
    using Mat4c   = Eigen::Matrix4cd;
    using Vec2c   = Eigen::Vector2cd;
    using Vec4c   = Eigen::Vector4cd;
    using Vec4d   = Eigen::Vector4d;

    inline void dirac_matrices_in(double xmn_in, bool cc_in);
    inline void define_spinors();
    void compute_hadron_tensor_eigen(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x, double q_y, double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4]);
    void compute_hadron_tensor_eigen_cc(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x, double q_y, double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4]);



    inline void sigccc(Mat4c &resp,
        double ff1v, double ff2v,
        double ffa,  double ffp);
    inline void shift(Mat4c &resp);

  }
}

#endif
