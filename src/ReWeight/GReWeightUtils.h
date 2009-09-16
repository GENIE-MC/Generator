
//____________________________________________________________________________
/*!

\namespace genie::utils::rew

\brief     Event reweighting utilities

\author    Jim Dobson <J.Dobson07 \at imperial.ac.uk>
           Imperial College London

           Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   Sep 09, 2009

\cpright   Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RW_UTILS_H_
#define _RW_UTILS_H_

#include "ReWeight/GSyst.h"

namespace genie {
namespace utils {
namespace rew   {

  //! Returns a weight to account for a change in hadron mean free path

  double MeanFreePathWeight(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double mfp_scale_factor, bool interacted,
    double nRpi=0.5, double nRnuc=1.0, double NR=3, double R0=1.4);

  double MeanFreePathWeight(
      double prob_def, double prob_twk, bool interacted);

  //!
  double FateXSec(genie::rew::GSyst_t syst, double kinE);


}  // rew   namespace
}  // utils namespace
}  // genie namespace

#endif // _RW_UTILS_H_
