//____________________________________________________________________________
/*!

\brief    Fractional uncertainties for physics parameters

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2008
*/
//____________________________________________________________________________

#ifndef _GSYST_ERROR_CONST_H_
#define _GSYST_ERROR_CONST_H_

namespace genie {
namespace rew {
  
  //
  // GENIE cross section parameter fractional uncertainties
  //
  const double kNuXSec_MaQEL_1SigmaErr        = 0.15;    ///< Ma QEL sigma 
  const double kNuXSec_MvQEL_1SigmaErr        = 0.05;    ///< Mv QEL sigma
  const double kNuXSec_MaRES_1SigmaErr        = 0.2;     ///< Ma RES sigma
  const double kNuXSec_MvRES_1SigmaErr        = 0.05;    ///< Mv RES sigma
  const double kNuXSec_MaCOHPi_1SigmaErr      = 0.1;     ///< Ma COHPi sigma
  const double kNuXSec_RvpCC1pi_1SigmaErr     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+p CC
  const double kNuXSec_RvpCC2pi_1SigmaErr     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+p CC
  const double kNuXSec_RvpNC1pi_1SigmaErr     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+p NC
  const double kNuXSec_RvpNC2pi_1SigmaErr     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+p NC
  const double kNuXSec_RvnCC1pi_1SigmaErr     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+n CC
  const double kNuXSec_RvnCC2pi_1SigmaErr     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+n CC
  const double kNuXSec_RvnNC1pi_1SigmaErr     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+n NC
  const double kNuXSec_RvnNC2pi_1SigmaErr     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+n NC
  const double kNuXSec_RvbarpCC1pi_1SigmaErr  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+p CC
  const double kNuXSec_RvbarpCC2pi_1SigmaErr  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+p CC
  const double kNuXSec_RvbarpNC1pi_1SigmaErr  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+p NC
  const double kNuXSec_RvbarpNC2pi_1SigmaErr  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+p NC
  const double kNuXSec_RvbarnCC1pi_1SigmaErr  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+n CC
  const double kNuXSec_RvbarnCC2pi_1SigmaErr  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+n CC
  const double kNuXSec_RvbarnNC1pi_1SigmaErr  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+n NC
  const double kNuXSec_RvbarnNC2pi_1SigmaErr  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+n NC

  //
  // GENIE INTRANUKE/hA parameter fractional uncertainties
  //

  // pi+,pi-,pi0
  const double kINuke_MFPTwk_pi_1SigmaErr       = 0.1;    ///< mean free path sigma
  const double kINuke_CExTwk_pi_1SigmaErr       = 0.1;    ///< charge exchange sigma
  const double kINuke_ElTwk_pi_1SigmaErr        = 0.1;    ///< elastic sigma
  const double kINuke_InelTwk_pi_1SigmaErr      = 0.1;    ///< inelastic sigma
  const double kINuke_AbsTwk_pi_1SigmaErr       = 0.1;    ///< absorption sigma
  const double kINuke_PiProdTwk_pi_1SigmaErr    = 0.1;    ///< pi production sigma

  // p, n
  const double kINuke_MFPTwk_N_1SigmaErr       = 0.1;    ///< mean free path sigma
  const double kINuke_CExTwk_N_1SigmaErr       = 0.1;    ///< charge exchange sigma
  const double kINuke_ElTwk_N_1SigmaErr        = 0.1;    ///< elastic sigma
  const double kINuke_InelTwk_N_1SigmaErr      = 0.1;    ///< inelastic sigma
  const double kINuke_AbsTwk_N_1SigmaErr       = 0.1;    ///< absorption sigma
  const double kINuke_PiProdTwk_N_1SigmaErr    = 0.1;    ///< pi production sigma

} // rew
} // genie

#endif


