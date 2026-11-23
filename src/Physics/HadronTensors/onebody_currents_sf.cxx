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
 *  \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
 *            For the full text of the license visit http://copyright.genie-mc.org
 *            or see $GENIE/LICENSE
 *
 * =====================================================================================
 */
#include <complex>
#include <array>
#include <cmath>
#include <iostream>

#include "Physics/HadronTensors/onebody_currents_sf.h"
#include "Physics/HadronTensors/TensorUtil.h"

namespace genie {
  namespace onebody_currents_sf {

    //===========================================================================
    // File-scope state (formerly Fortran SAVE variables) and helpers.
    //===========================================================================
    namespace {

      // Constants -----------------------------------------------------------------
      const std::complex<double> kCZero{0.0, 0.0};
      const std::complex<double> kCOne {1.0, 0.0};
      const std::complex<double> kCI   {0.0, 1.0};

      // Pauli matrices, identities, Dirac matrices, etc. -------------------------
      genie::TensorUtil::Matrix2cd gSigma[3];
      genie::TensorUtil::Matrix2cd gId2;
      genie::TensorUtil::Matrix4cd gId4;

      genie::TensorUtil::Matrix4cd gGammaMu[5];          // gamma^0..gamma^3, gamma^5
      genie::TensorUtil::Matrix4cd gMetric;              // metric diag(+1,-1,-1,-1)
      genie::TensorUtil::Matrix4cd gSigmaMuNu[4][4];     // sigma^{mu nu}, mu, nu = 0..3

      genie::TensorUtil::Matrix4cd gQSlash;              // unused but kept for completeness

      // Spinors (2-component and 4-component) ------------------------------------
      genie::TensorUtil::Vector2cd gUp;
      genie::TensorUtil::Vector2cd gDown;
      std::array<genie::TensorUtil::Vector4cd, 2> gUp1;
      std::array<genie::TensorUtil::Vector4cd, 2> gUpp1;
      std::array<genie::TensorUtil::Vector4cd, 2> gUbarP1;
      std::array<genie::TensorUtil::Vector4cd, 2> gUbarPp1;

      // 4-vectors and scalar -----------------------------------------------------
      genie::TensorUtil::Vector4d gP1;
      genie::TensorUtil::Vector4d gPp1;
      genie::TensorUtil::Vector4d gQt;
      double          gW = 0.0;

      // Current tensor J_1^{mu} (4x4 matrix for each mu = 0..3) ------------------
      genie::TensorUtil::Matrix4cd gJ1[4];

      // Nucleon mass and current-conservation flag -------------------------------
      double gXmn = 0.0;
      bool   gCC  = false;


      //---------------------------------------------------------------------------
      // DiracMatricesIn: initialize gamma, sigma, metric, etc.
      //
      // Fortran:
      //   subroutine dirac_matrices_in(xmn_in, cc_in)
      //---------------------------------------------------------------------------
      void DiracMatricesIn(double xmn_in, bool cc_in)
      {
        gXmn = xmn_in;
        gCC  = cc_in;

        // 2x2 identity
        gId2.setZero();
        gId2(0,0) = kCOne;
        gId2(1,1) = kCOne;

        // Pauli matrices
        gSigma[0].setZero();
        gSigma[0](0,1) = kCOne;
        gSigma[0](1,0) = kCOne;

        gSigma[1].setZero();
        gSigma[1](0,1) = -kCI;
        gSigma[1](1,0) =  kCI;

        gSigma[2].setZero();
        gSigma[2](0,0) =  kCOne;
        gSigma[2](1,1) = -kCOne;

        // gamma matrices
        for(int k = 0; k < 5; ++k)
        {
          gGammaMu[k].setZero();
        }

        // gamma^0: diag(I, -I)
        gGammaMu[0].block<2,2>(0,0) =  gId2;
        gGammaMu[0].block<2,2>(2,2) = -gId2;

        // id4: diag(I2, I2)
        gId4.setZero();
        gId4.block<2,2>(0,0) = gId2;
        gId4.block<2,2>(2,2) = gId2;

        // gamma^i for i = 1..3
        for(int mu = 1; mu <= 3; ++mu)
        {
          int s = mu - 1;
          gGammaMu[mu].block<2,2>(0,2) =  gSigma[s];
          gGammaMu[mu].block<2,2>(2,0) = -gSigma[s];
        }

        // gamma^5
        gGammaMu[4].block<2,2>(0,2) = gId2;
        gGammaMu[4].block<2,2>(2,0) = gId2;

        // metric g^{mu nu} = diag(+1,-1,-1,-1)
        gMetric.setZero();
        gMetric(0,0) =  kCOne;
        gMetric(1,1) = -kCOne;
        gMetric(2,2) = -kCOne;
        gMetric(3,3) = -kCOne;

        // sigma^{mu nu} = i/2 [gamma^mu, gamma^nu]
        for(int mu = 0; mu < 4; ++mu)
        {
          for(int nu = 0; nu < 4; ++nu)
          {
            gSigmaMuNu[mu][nu] =
              (kCI * 0.5) * (gGammaMu[mu] * gGammaMu[nu] -
                  gGammaMu[nu] * gGammaMu[mu]);
          }
        }

        // 2-component spinors
        gUp   << kCOne, kCZero;
        gDown << kCZero, kCOne;

        // clear q_sl just in case
        gQSlash.setZero();
      }


      //---------------------------------------------------------------------------
      // DefineSpinors: build 4-component spinors from gP1, gPp1, gXmn.
      //
      // Fortran:
      //   subroutine define_spinors()
      //---------------------------------------------------------------------------
      void DefineSpinors()
      {
        genie::TensorUtil::Matrix2cd sigp1  = genie::TensorUtil::Matrix2cd::Zero();
        genie::TensorUtil::Matrix2cd sigpp1 = genie::TensorUtil::Matrix2cd::Zero();

        for(int s = 0; s < 2; ++s)
        {
          gUp1     [s].setZero();
          gUpp1    [s].setZero();
          gUbarP1  [s].setZero();
          gUbarPp1 [s].setZero();
        }

        double cp1  = std::sqrt((gP1 (0) + gXmn) / 2.0);
        double cpp1 = std::sqrt((gPp1(0) + gXmn) / 2.0);

        for(int i = 0; i < 3; ++i)
        {
          sigp1  += gSigma[i] * gP1 (i + 1);
          sigpp1 += gSigma[i] * gPp1(i + 1);
        }

        // ---- gUp1 (incoming quadrispinor) --------------------------------------
        {
          double denom = gP1(0) + gXmn;

          gUp1[0].head<2>() = gUp;
          genie::TensorUtil::Vector2cd tmp = sigp1 * gUp;
          gUp1[0].segment<2>(2) = tmp / denom;

          gUp1[1].head<2>() = gDown;
          tmp = sigp1 * gDown;
          gUp1[1].segment<2>(2) = tmp / denom;

          gUp1[0] *= cp1;
          gUp1[1] *= cp1;
        }

        // ---- gUpp1 (outgoing quadrispinor) -------------------------------------
        {
          double denom = gPp1(0) + gXmn;

          gUpp1[0].head<2>() = gUp;
          genie::TensorUtil::Vector2cd tmp = sigpp1 * gUp;
          gUpp1[0].segment<2>(2) = tmp / denom;

          gUpp1[1].head<2>() = gDown;
          tmp = sigpp1 * gDown;
          gUpp1[1].segment<2>(2) = tmp / denom;

          gUpp1[0] *= cpp1;
          gUpp1[1] *= cpp1;
        }

        // ---- gUbarP1 (incoming adjoint) ----------------------------------------
        {
          double denom = gP1(0) + gXmn;

          gUbarP1[0].head<2>() = gUp;
          genie::TensorUtil::RowVector2cd up_row  = gUp.transpose();
          genie::TensorUtil::RowVector2cd tmp_row = up_row * sigp1;
          gUbarP1[0](2) = -tmp_row(0) / denom;
          gUbarP1[0](3) = -tmp_row(1) / denom;

          gUbarP1[1].head<2>() = gDown;
          genie::TensorUtil::RowVector2cd down_row = gDown.transpose();
          tmp_row = down_row * sigp1;
          gUbarP1[1](2) = -tmp_row(0) / denom;
          gUbarP1[1](3) = -tmp_row(1) / denom;

          gUbarP1[0] *= cp1;
          gUbarP1[1] *= cp1;
        }

        // ---- gUbarPp1 (outgoing adjoint) ---------------------------------------
        {
          double denom = gPp1(0) + gXmn;

          gUbarPp1[0].head<2>() = gUp;
          genie::TensorUtil::RowVector2cd up_row  = gUp.transpose();
          genie::TensorUtil::RowVector2cd tmp_row = up_row * sigpp1;
          gUbarPp1[0](2) = -tmp_row(0) / denom;
          gUbarPp1[0](3) = -tmp_row(1) / denom;

          gUbarPp1[1].head<2>() = gDown;
          genie::TensorUtil::RowVector2cd down_row = gDown.transpose();
          tmp_row = down_row * sigpp1;
          gUbarPp1[1](2) = -tmp_row(0) / denom;
          gUbarPp1[1](3) = -tmp_row(1) / denom;

          gUbarPp1[0] *= cpp1;
          gUbarPp1[1] *= cpp1;
        }
      }


      //---------------------------------------------------------------------------
      // CurrentInit: copy 4-momenta and w into module state.
      //
      // Fortran:
      //   subroutine current_init(p1_in, pp1_in, qt_in, w_in)
      //---------------------------------------------------------------------------
      void CurrentInit(const genie::TensorUtil::Vector4d & p1_in,
          const genie::TensorUtil::Vector4d & pp1_in,
          const genie::TensorUtil::Vector4d & qt_in,
          double w_in)
      {
        gP1  = p1_in;
        gPp1 = pp1_in;
        gQt  = qt_in;
        gW   = w_in;
      }


      //---------------------------------------------------------------------------
      // DetJa: build J_1(mu) = J_V(mu) + J_A(mu).
      //
      // Fortran:
      //   subroutine det_Ja(f1v, f2v, ffa, ffp)
      //---------------------------------------------------------------------------
      void DetJa(double f1v, double f2v, double ffa, double ffp)
      {
        genie::TensorUtil::Matrix4cd J_1_V[4];
        genie::TensorUtil::Matrix4cd J_1_A[4];

        for(int mu = 0; mu < 4; ++mu)
        {
          J_1_V[mu].setZero();
          J_1_A[mu].setZero();

          // Vector part: sum over nu
          for(int nu = 0; nu < 4; ++nu)
          {
            std::complex<double> g      = gMetric(nu, nu);
            double               q      = gQt(nu);
            std::complex<double> factor = kCI * f2v * g * q / (2.0 * gXmn);
            J_1_V[mu] += factor * gSigmaMuNu[mu][nu];
          }

          J_1_V[mu] += f1v * gGammaMu[mu];

          // Axial part: ffa * gamma_mu(mu) * gamma5
          genie::TensorUtil::Matrix4cd tmp = gGammaMu[mu] * gGammaMu[4];
          J_1_A[mu] += ffa * tmp;

          // Pseudoscalar piece: ffp * gamma5 * q_mu / xmn
          J_1_A[mu] += ffp * (gQt(mu) / gXmn) * gGammaMu[4];
        }

        if(gCC)
        {
          double factor = gW / gQt(3);
          J_1_V[3] = factor * J_1_V[0];
        }

        for(int mu = 0; mu < 4; ++mu)
        {
          gJ1[mu] = J_1_V[mu] + J_1_A[mu];
        }
      }


      //---------------------------------------------------------------------------
      // DetRes1b: build hadronic tensor from gJ1 and spinors.
      //
      // Fortran:
      //   subroutine det_res1b(res)
      //---------------------------------------------------------------------------
      void DetRes1b(genie::TensorUtil::Matrix4cd & res)
      {
        std::array<std::array<genie::TensorUtil::Vector4cd, 2>, 2> J_mu;

        for(int i1 = 0; i1 < 2; ++i1)
        {
          for(int f1 = 0; f1 < 2; ++f1)
          {
            for(int mu = 0; mu < 4; ++mu)
            {
              genie::TensorUtil::Vector4cd tmp = gJ1[mu] * gUp1[i1];

              // NOTE: Fortran sum(a*b) has NO implicit conjugation.
              std::complex<double> val = kCZero;
              for(int a = 0; a < 4; ++a)
              {
                val += gUbarPp1[f1](a) * tmp(a);
              }

              J_mu[f1][i1](mu) = val;
            }
          }
        }

        // res(i,j) = sum_{spins, mu} J_mu^*(mu)_i * J_mu(mu)_j
        res.setZero();
        for(int i1 = 0; i1 < 2; ++i1)
        {
          for(int f1 = 0; f1 < 2; ++f1)
          {
            const genie::TensorUtil::Vector4cd & Jv = J_mu[f1][i1];
            for(int i = 0; i < 4; ++i)
            {
              for(int j = 0; j < 4; ++j)
              {
                res(i,j) += std::conj(Jv(i)) * Jv(j);
              }
            }
          }
        }
      }


      //---------------------------------------------------------------------------
      // Shift: CSHIFT(resp,1,DIM=1) then CSHIFT(resp,1,DIM=2).
      //
      // Fortran: subroutine shift(resp) bind(C, name="shift2")
      //---------------------------------------------------------------------------
      void Shift(genie::TensorUtil::Matrix4cd & resp)
      {
        genie::TensorUtil::Matrix4cd tmp;

        for(int i = 0; i < 4; ++i)
        {
          int src_i = (i + 1 + 4) % 4;
          for(int j = 0; j < 4; ++j)
          {
            tmp(i, j) = resp(src_i, j);
          }
        }
        resp = tmp;

        for(int j = 0; j < 4; ++j)
        {
          int src_j = (j + 1 + 4) % 4;
          for(int i = 0; i < 4; ++i)
          {
            tmp(i, j) = resp(i, src_j);
          }
        }
        resp = tmp;
      }


      //---------------------------------------------------------------------------
      // SigCCC: DetJa -> DetRes1b -> spin average -> Shift -> transpose.
      //
      // Fortran: subroutine sigccc(resp,...) bind(C, name="sigccc2")
      //---------------------------------------------------------------------------
      void SigCCC(genie::TensorUtil::Matrix4cd & resp,
          double ff1v, double ff2v,
          double ffa,  double ffp)
      {
        DetJa(ff1v, ff2v, ffa, ffp);

        genie::TensorUtil::Matrix4cd tmp;
        DetRes1b(tmp);

        resp = 0.5 * tmp;
        Shift(resp);
        resp.transposeInPlace();
      }


      //---------------------------------------------------------------------------
      // Common implementation for both ComputeHadronTensor entry points.
      //---------------------------------------------------------------------------
      void ComputeHadronTensorImpl(double xmn_in, double w_in, double wt,
          double xk_x, double xk_y, double xk_z,
          double q_x,  double q_y,  double q_z,
          double ff1v, double ff2v,
          double ffa,  double ffp,
          bool conserve_current,
          std::complex<double> HadronTensor[4][4])
      {
        genie::TensorUtil::Matrix4cd resp;
        resp.setZero();

        // Set up Dirac matrices and module globals (xmn, cc, etc.)
        DiracMatricesIn(xmn_in, conserve_current);

        // |k| and |k+q|
        double xk = std::sqrt(xk_x*xk_x + xk_y*xk_y + xk_z*xk_z);
        double xp = std::sqrt((xk_x + q_x)*(xk_x + q_x) +
            (xk_y + q_y)*(xk_y + q_y) +
            (xk_z + q_z)*(xk_z + q_z));

        // Energies: Fortran uses module variable xmn (not xmn_in) here
        double ek  = std::sqrt(gXmn*gXmn + xk*xk);
        double epf = std::sqrt(gXmn*gXmn + xp*xp);

        genie::TensorUtil::Vector4d qt_4, p_4, pp_4;
        qt_4 << wt, q_x, q_y, q_z;
        p_4  << ek, xk_x, xk_y, xk_z;

        pp_4(0)        = epf;
        pp_4.tail<3>() = p_4.tail<3>() + qt_4.tail<3>();

        CurrentInit(p_4, pp_4, qt_4, w_in);
        DefineSpinors();

        SigCCC(resp, ff1v, ff2v, ffa, ffp);

        for(int i = 0; i < 4; ++i)
        {
          for(int j = 0; j < 4; ++j)
          {
            HadronTensor[i][j] = resp(j, i);
          }
        }
      }

    }  // anonymous namespace


    //---------------------------------------------------------------------------
    // ComputeNucleonTensor: entry point WITHOUT current conservation.
    //---------------------------------------------------------------------------
    void ComputeNucleonTensor(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x,  double q_y,  double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4])
    {
      // conserve_current = false (matches original behavior)
      ComputeHadronTensorImpl(xmn_in, w_in, wt,
          xk_x, xk_y, xk_z,
          q_x,  q_y,  q_z,
          ff1v, ff2v, ffa, ffp,
          false,
          HadronTensor);
    }


    //---------------------------------------------------------------------------
    // ComputeNucleonTensorCC: entry point WITH current conservation.
    //---------------------------------------------------------------------------
    void ComputeNucleonTensorCC(double xmn_in, double w_in, double wt,
        double xk_x, double xk_y, double xk_z,
        double q_x,  double q_y,  double q_z,
        double ff1v, double ff2v,
        double ffa,  double ffp,
        std::complex<double> HadronTensor[4][4])
    {
      // conserve_current = true
      ComputeHadronTensorImpl(xmn_in, w_in, wt,
          xk_x, xk_y, xk_z,
          q_x,  q_y,  q_z,
          ff1v, ff2v, ffa, ffp,
          true,
          HadronTensor);
    }

  }  // namespace onebody_currents_sf
}  // namespace genie
