//____________________________________________________________________________
/*!

\brief    Fractional uncertainties for physics parameters

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_FRACERRORS_H_
#define _G_FRACERRORS_H_

namespace genie {
namespace rew   {

  //
  // GENIE cross section parameter fractional uncertainties
  //

  const double kNuXSec_MaQELFracError        = 0.15;    ///< Ma QEL sigma 
  const double kNuXSec_MvQELFracError        = 0.05;    ///< Mv QEL sigma
  const double kNuXSec_MaRESFracError        = 0.2;     ///< Ma RES sigma
  const double kNuXSec_MvRESFracError        = 0.05;    ///< Mv RES sigma
  const double kNuXSec_MaCOHPiError          = 0.1;     ///< Ma COHPi sigma
  const double kNuXSec_RvpCC1piFracError     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+p CC
  const double kNuXSec_RvpCC2piFracError     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+p CC
  const double kNuXSec_RvpNC1piFracError     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+p NC
  const double kNuXSec_RvpNC2piFracError     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+p NC
  const double kNuXSec_RvnCC1piFracError     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+n CC
  const double kNuXSec_RvnCC2piFracError     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+n CC
  const double kNuXSec_RvnNC1piFracError     = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for v+n NC
  const double kNuXSec_RvnNC2piFracError     = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for v+n NC
  const double kNuXSec_RvbarpCC1piFracError  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+p CC
  const double kNuXSec_RvbarpCC2piFracError  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+p CC
  const double kNuXSec_RvbarpNC1piFracError  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+p NC
  const double kNuXSec_RvbarpNC2piFracError  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+p NC
  const double kNuXSec_RvbarnCC1piFracError  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+n CC
  const double kNuXSec_RvbarnCC2piFracError  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+n CC
  const double kNuXSec_RvbarnNC1piFracError  = 0.5;     ///< sigma for 1pi non-RES bkg in the RES region, for vbar+n NC
  const double kNuXSec_RvbarnNC2piFracError  = 0.5;     ///< sigma for 2pi non-RES bkg in the RES region, for vbar+n NC

  
  //
  // GENIE INTRANUKE/hA parameter fractional uncertainties
  //

  //
  // These are simplest possible errors, fixed %, which will be replaced
  // with hadron energy dependent splines in future.
  //

  // pi+,pi-,pi0
  const double kPionMFPFracError             = 0.1;    ///< mean free path sigma
  const double kPion_CExFracError            = 0.1;    ///< charge exchange sigma
  const double kPion_ElFracError             = 0.1;    ///< elastic sigma
  const double kPion_InelFracError           = 0.1;    ///< inelastic sigma
  const double kPion_AbsFracError            = 0.1;    ///< absorption sigma
  const double kPion_PiProdFracError         = 0.1;    ///< pi production sigma

  // p, n
  const double kNuclMFPFracError             = 0.1;    ///< mean free path sigma
  const double kNucl_CExFracError            = 0.1;    ///< charge exchange sigma
  const double kNucl_ElFracError             = 0.1;    ///< elastic sigma
  const double kNucl_InelFracError           = 0.1;    ///< inelastic sigma
  const double kNucl_AbsFracError            = 0.1;    ///< absorption sigma
  const double kNucl_PiProdFracError         = 0.1;    ///< pi production sigma

}  // rew   namespace    
}  // genie namespace

#endif


