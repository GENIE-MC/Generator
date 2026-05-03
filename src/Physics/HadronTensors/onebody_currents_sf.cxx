/*
 * =====================================================================================
 *
 *       Filename:  onebody_currents_sf.cxx
 *
 *    Description: Native C++ implementation, replacing the previous Fortran wrapper.
 *
 *                 Reference: https://arxiv.org/pdf/2308.15524
 *
 *        Version:  1.0
 *        Created:  12/09/2025 10:45:26 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */

#include <Eigen/Dense>
#include <complex>
#include <array>
#include <cmath>
#include <iostream>
#include "Physics/HadronTensors/onebody_currents_sf.h"
namespace genie{
  namespace onebody_currents_sf {

    // ----------------------------------------------------------------------
    // "Module" variables (Fortran SAVE variables)
    // ----------------------------------------------------------------------

    // constants
    inline constexpr cdouble czero{0.0, 0.0};
    inline constexpr cdouble cone {1.0, 0.0};
    inline constexpr cdouble ci   {0.0, 1.0};

    // Pauli matrices, identity, Dirac matrices, etc.
    inline Mat2c sigma[3];          // 3 Pauli matrices
    inline Mat2c id2;               // 2x2 identity
    inline Mat4c id4;               // 4x4 identity

    inline Mat4c gamma_mu[5];       // gamma^0..gamma^3, gamma^5
    inline Mat4c g_munu;            // metric diag(+1,-1,-1,-1)
    inline Mat4c sigma_munu[4][4];  // sigma^{mu nu}, mu,nu = 0..3

    inline Mat4c q_sl;              // unused but kept for completeness

    // spinors (2-component and 4-component)
    inline Vec2c up, down;
    inline std::array<Vec4c, 2> up1;      // up1(spin,dirac)
    inline std::array<Vec4c, 2> upp1;     // upp1(spin,dirac)
    inline std::array<Vec4c, 2> ubarp1;   // ubarp1(spin,dirac)
    inline std::array<Vec4c, 2> ubarpp1;  // ubarpp1(spin,dirac)

    // 4-vectors and scalar
    inline Vec4d p1, pp1, qt;
    inline double w = 0.0;

    // current tensor J_1^{mu} (4x4 matrix for each mu=0..3)
    inline Mat4c J_1[4];

    // nucleon mass (public in Fortran)
    inline double xmn = 0.0;

    // current conservation flag
    inline bool cc = false;

    // ----------------------------------------------------------------------
    // dirac_matrices_in: initialize gamma, sigma, metric, etc.
    // ----------------------------------------------------------------------
    //
    // Fortran:
    //   subroutine dirac_matrices_in(xmn_in, cc_in)
    // ----------------------------------------------------------------------
    inline void dirac_matrices_in(double xmn_in, bool cc_in)
    {
      xmn = xmn_in;
      cc  = cc_in;

      // 2x2 identity
      id2.setZero();
      id2(0,0) = cone;
      id2(1,1) = cone;

      // Pauli matrices (Fortran sig(1..3,:,:))
      // sigma(1)
      sigma[0].setZero();
      sigma[0](0,1) = cone;
      sigma[0](1,0) = cone;
      // sigma(2)
      sigma[1].setZero();
      sigma[1](0,1) = -ci;
      sigma[1](1,0) =  ci;
      // sigma(3)
      sigma[2].setZero();
      sigma[2](0,0) =  cone;
      sigma[2](1,1) = -cone;

      // gamma matrices
      for (int k = 0; k < 5; ++k)
        gamma_mu[k].setZero();

      // gamma^0: diag(I, -I)
      gamma_mu[0].setZero();
      gamma_mu[0].block<2,2>(0,0) = id2;
      gamma_mu[0].block<2,2>(2,2) = -id2;

      // id4: diag(I2, I2)
      id4.setZero();
      id4.block<2,2>(0,0) = id2;
      id4.block<2,2>(2,2) = id2;

      // gamma^i for i = 1..3
      for (int mu = 1; mu <= 3; ++mu) {
        gamma_mu[mu].setZero();
        int s = mu - 1; // sigma index 0..2
        gamma_mu[mu].block<2,2>(0,2) =  sigma[s];
        gamma_mu[mu].block<2,2>(2,0) = -sigma[s];
      }

      // gamma^5
      gamma_mu[4].setZero();
      gamma_mu[4].block<2,2>(0,2) = id2;
      gamma_mu[4].block<2,2>(2,0) = id2;

      // metric g^{mu nu} = diag(+1,-1,-1,-1)
      g_munu.setZero();
      g_munu(0,0) = cone;
      g_munu(1,1) = -cone;
      g_munu(2,2) = -cone;
      g_munu(3,3) = -cone;

      // sigma^{mu nu} = i/2 [gamma^mu, gamma^nu], mu,nu = 0..3
      for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
          sigma_munu[mu][nu] =
            (ci * 0.5) * (gamma_mu[mu] * gamma_mu[nu] -
                gamma_mu[nu] * gamma_mu[mu]);
        }
      }

      // 2-component spinors
      up   << cone, czero;
      down << czero, cone;

      // clear q_sl just in case
      q_sl.setZero();
    }

    // ----------------------------------------------------------------------
    // define_spinors: build 4-component spinors from p1, pp1, xmn
    // ----------------------------------------------------------------------
    //
    // Fortran:
    //   subroutine define_spinors()
    // ----------------------------------------------------------------------
    inline void define_spinors()
    {
      Mat2c sigp1  = Mat2c::Zero();
      Mat2c sigpp1 = Mat2c::Zero();

      // zero spinor arrays
      for (int s = 0; s < 2; ++s) {
        up1[s].setZero();
        upp1[s].setZero();
        ubarp1[s].setZero();
        ubarpp1[s].setZero();
      }

      // normalization factors
      double cp1  = std::sqrt((p1(0)  + xmn) / 2.0);
      double cpp1 = std::sqrt((pp1(0) + xmn) / 2.0);

      // sigp1 = sum_i sigma_i * p1(i+1), i=0..2 → p1(1..3)
      // sigpp1 similarly with pp1
      for (int i = 0; i < 3; ++i) {
        sigp1  += sigma[i] * p1(i + 1);
        sigpp1 += sigma[i] * pp1(i + 1);
      }

      // ---- build up1 (incoming quadrispinor) ---------------------------
      {
        double denom = p1(0) + xmn;

        // spin up
        up1[0].head<2>() = up;                     // components 0..1
        Vec2c tmp = sigp1 * up;
        up1[0].segment<2>(2) = tmp / denom;        // components 2..3

        // spin down
        up1[1].head<2>() = down;
        tmp = sigp1 * down;
        up1[1].segment<2>(2) = tmp / denom;

        // overall normalization
        up1[0] *= cp1;
        up1[1] *= cp1;
      }

      // ---- build upp1 (outgoing quadrispinor) -------------------------
      {
        double denom = pp1(0) + xmn;

        // spin up
        upp1[0].head<2>() = up;
        Vec2c tmp = sigpp1 * up;
        upp1[0].segment<2>(2) = tmp / denom;

        // spin down
        upp1[1].head<2>() = down;
        tmp = sigpp1 * down;
        upp1[1].segment<2>(2) = tmp / denom;

        // overall normalization
        upp1[0] *= cpp1;
        upp1[1] *= cpp1;
      }

      // ---- build ubarp1 (incoming adjoint) ----------------------------
      {
        double denom = p1(0) + xmn;

        // spin up
        ubarp1[0].head<2>() = up;
        Eigen::RowVector2cd up_row  = up.transpose();
        Eigen::RowVector2cd tmp_row = up_row * sigp1;  // row * matrix
        ubarp1[0](2) = -tmp_row(0) / denom;
        ubarp1[0](3) = -tmp_row(1) / denom;

        // spin down
        ubarp1[1].head<2>() = down;
        Eigen::RowVector2cd down_row = down.transpose();
        tmp_row = down_row * sigp1;
        ubarp1[1](2) = -tmp_row(0) / denom;
        ubarp1[1](3) = -tmp_row(1) / denom;

        // overall normalization
        ubarp1[0] *= cp1;
        ubarp1[1] *= cp1;
      }

      // ---- build ubarpp1 (outgoing adjoint) ---------------------------
      {
        double denom = pp1(0) + xmn;

        // spin up
        ubarpp1[0].head<2>() = up;
        Eigen::RowVector2cd up_row  = up.transpose();
        Eigen::RowVector2cd tmp_row = up_row * sigpp1;
        ubarpp1[0](2) = -tmp_row(0) / denom;
        ubarpp1[0](3) = -tmp_row(1) / denom;

        // spin down
        ubarpp1[1].head<2>() = down;
        Eigen::RowVector2cd down_row = down.transpose();
        tmp_row = down_row * sigpp1;
        ubarpp1[1](2) = -tmp_row(0) / denom;
        ubarpp1[1](3) = -tmp_row(1) / denom;

        // overall normalization
        ubarpp1[0] *= cpp1;
        ubarpp1[1] *= cpp1;
      }
    }
    // ----------------------------------------------------------------------
    // current_init: copy 4-momenta and w into module state
    // ----------------------------------------------------------------------
    //
    // Fortran:
    //   subroutine current_init(p1_in, pp1_in, qt_in, w_in)
    // ----------------------------------------------------------------------
    inline void current_init(const Vec4d &p1_in,
        const Vec4d &pp1_in,
        const Vec4d &qt_in,
        double w_in)
    {
      p1  = p1_in;
      pp1 = pp1_in;
      qt  = qt_in;
      w   = w_in;
    }

    // ----------------------------------------------------------------------
    // det_Ja: build J_1(mu) = J_V(mu) + J_A(mu)
    // ----------------------------------------------------------------------
    //
    // Fortran:
    //   subroutine det_Ja(f1v, f2v, ffa, ffp)
    // ----------------------------------------------------------------------
    inline void det_Ja(double f1v, double f2v, double ffa, double ffp)
    {
      Mat4c J_1_V[4];
      Mat4c J_1_A[4];

      for (int mu = 0; mu < 4; ++mu) {
        J_1_V[mu].setZero();
        J_1_A[mu].setZero();

        // Vector part: sum over nu
        for (int nu = 0; nu < 4; ++nu) {
          cdouble g = g_munu(nu,nu); // +1 or -1
          double q = qt(nu);
          cdouble factor = ci * f2v * g * q / (2.0 * xmn);
          J_1_V[mu] += factor * sigma_munu[mu][nu];
        }

        // Add f1v * gamma_mu
        J_1_V[mu] += f1v * gamma_mu[mu];

        // Axial part: ffa * gamma_mu(mu) * gamma5
        Mat4c tmp = gamma_mu[mu] * gamma_mu[4]; // gamma^mu * gamma^5
        J_1_A[mu] += ffa * tmp;

        // Pseudoscalar piece: ffp * gamma5 * q_mu / xmn
        J_1_A[mu] += ffp * (qt(mu) / xmn) * gamma_mu[4];
      }

      // Current conservation (q0 J0 = q3 J3) when q along z
      if (cc) {
        // Fortran: J_1_V(:,:,4) = (w/qt(4))*J_1_V(:,:,1)
        // Indices: mu=3→gamma_mu[3], mu=0→gamma_mu[0], qt(4)→qt(3)
        double factor = w / qt(3);
        J_1_V[3] = factor * J_1_V[0];
      }

      // Store total currents J_1(mu)
      for (int mu = 0; mu < 4; ++mu)
        J_1[mu] = J_1_V[mu] + J_1_A[mu];
    }

    // ----------------------------------------------------------------------
    // det_res1b: build hadronic tensor from J_1 and spinors
    // ----------------------------------------------------------------------
    //
    // Fortran:
    //   subroutine det_res1b(res)
    // ----------------------------------------------------------------------
    inline void det_res1b(Mat4c &res)
    {
      // J_mu(f1, i1, mu) as in Fortran: spin_out, spin_in, Lorentz component
      std::array<std::array<Vec4c, 2>, 2> J_mu;

      // Fill J_mu
      for (int i1 = 0; i1 < 2; ++i1) {      // initial spin
        for (int f1 = 0; f1 < 2; ++f1) {  // final spin
          for (int mu = 0; mu < 4; ++mu) {
            // tmp = J_1(:,:,mu) * up1(i1,:)
            Vec4c tmp = J_1[mu] * up1[i1];

            // J_mu(f1,i1,mu) = sum_alpha ubarpp1(f1,alpha)*tmp(alpha)
            // NOTE: Fortran sum(a*b) has NO implicit conjugation.
            cdouble val = czero;
            for (int a = 0; a < 4; ++a)
              val += ubarpp1[f1](a) * tmp(a);

            J_mu[f1][i1](mu) = val;
          }
        }
      }

      // res(i,j) = sum_{spins,mu} J_mu^†(mu)_i * J_mu(mu)_j
      res.setZero();
      for (int i1 = 0; i1 < 2; ++i1) {
        for (int f1 = 0; f1 < 2; ++f1) {
          const Vec4c &Jv = J_mu[f1][i1];
          for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
              res(i,j) += std::conj(Jv(i)) * Jv(j);
            }
          }
        }
      }
    }

    void compute_hadron_tensor_eigen(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x, double q_y, double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4]){

      Eigen::Matrix4cd resp; resp.setZero();
      // Use current conservation (conserve_current = .TRUE.)
      bool conserve_current = false;

      // Set up Dirac matrices and module globals (xmn, cc, etc.)
      dirac_matrices_in(xmn_in, conserve_current);

      // |k| and |k+q|
      double xk = std::sqrt(xk_x*xk_x + xk_y*xk_y + xk_z*xk_z);
      double xp = std::sqrt((xk_x+q_x)*(xk_x+q_x) +
          (xk_y+q_y)*(xk_y+q_y) +
          (xk_z+q_z)*(xk_z+q_z));

      // Energies: Fortran uses module variable xmn (not xmn_in) here
      double ek  = std::sqrt(xmn * xmn + xk * xk);
      double epf = std::sqrt(xmn * xmn + xp * xp);

      // 4-vectors: (t,x,y,z) = (E,px,py,pz) and (wt,qx,qy,qz)
      Vec4d qt_4, p_4, pp_4;
      qt_4 << wt, q_x, q_y, q_z;

      p_4 << ek, xk_x, xk_y, xk_z;

      pp_4(0) = epf;
      pp_4.tail<3>() = p_4.tail<3>() + qt_4.tail<3>();

      // Initialize current and spinors
      current_init(p_4, pp_4, qt_4, w_in);
      define_spinors();

      // Build response tensor
      sigccc(resp, ff1v, ff2v, ffa, ffp);

      // fill the C array
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          HadronTensor[i][j] = resp(j, i);
        }
      }

    }


    // current conservation
    void compute_hadron_tensor_eigen_cc(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x, double q_y, double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4]){

      Eigen::Matrix4cd resp; resp.setZero();
      // Use current conservation (conserve_current = .TRUE.)
      bool conserve_current = true;

      // Set up Dirac matrices and module globals (xmn, cc, etc.)
      dirac_matrices_in(xmn_in, conserve_current);

      // |k| and |k+q|
      double xk = std::sqrt(xk_x*xk_x + xk_y*xk_y + xk_z*xk_z);
      double xp = std::sqrt((xk_x+q_x)*(xk_x+q_x) +
          (xk_y+q_y)*(xk_y+q_y) +
          (xk_z+q_z)*(xk_z+q_z));

      // Energies: Fortran uses module variable xmn (not xmn_in) here
      double ek  = std::sqrt(xmn * xmn + xk * xk);
      double epf = std::sqrt(xmn * xmn + xp * xp);

      // 4-vectors: (t,x,y,z) = (E,px,py,pz) and (wt,qx,qy,qz)
      Vec4d qt_4, p_4, pp_4;
      qt_4 << wt, q_x, q_y, q_z;

      p_4 << ek, xk_x, xk_y, xk_z;

      pp_4(0) = epf;
      pp_4.tail<3>() = p_4.tail<3>() + qt_4.tail<3>();

      // Initialize current and spinors
      current_init(p_4, pp_4, qt_4, w_in);
      define_spinors();

      // Build response tensor
      sigccc(resp, ff1v, ff2v, ffa, ffp);

      // fill the C array
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          HadronTensor[i][j] = resp(j, i);
        }
      }

    }




    // ----------------------------------------------------------------------
    // shift2: CSHIFT(resp,1,DIM=1) then CSHIFT(resp,1,DIM=2)
    // Fortran: subroutine shift(resp) bind(C,name="shift2")
    // ----------------------------------------------------------------------
    inline void shift(Mat4c &resp)
    {
      Mat4c tmp;

      // CSHIFT along DIM=1 (first index, i -> i+1 with wrap)
      for (int i = 0; i < 4; ++i) {
        int src_i = (i + 1 + 4) % 4;  // old(i-1, j)
        for (int j = 0; j < 4; ++j) {
          tmp(i,j) = resp(src_i,j);
        }
      }
      resp = tmp;

      // CSHIFT along DIM=2 (second index, j -> j+1 with wrap)
      for (int j = 0; j < 4; ++j) {
        int src_j = (j + 1 + 4) % 4;  // old(i, j-1)
        for (int i = 0; i < 4; ++i) {
          tmp(i,j) = resp(i,src_j);
        }
      }
      resp = tmp;
    }

    // ----------------------------------------------------------------------
    // sigccc2: det_Ja -> det_res1b -> spin average -> shift2 -> transpose
    // Fortran: subroutine sigccc(resp,...) bind(C,name="sigccc2")
    // ----------------------------------------------------------------------
    inline void sigccc(Mat4c &resp,
        double ff1v, double ff2v,
        double ffa,  double ffp)
    {
      // Build current J_1
      det_Ja(ff1v, ff2v, ffa, ffp);

      // Compute response tensor
      Mat4c tmp;
      det_res1b(tmp);

      // Spin average (factor 1/2)
      resp = 0.5 * tmp;

      // Shift (t,x,y,z) -> (x,y,z,t)
      shift(resp);

      // Transpose: inter = TRANSPOSE(resp); resp = inter
      resp.transposeInPlace();

    }
  } // namespace onebody_currents_sf
} // namespace genie

