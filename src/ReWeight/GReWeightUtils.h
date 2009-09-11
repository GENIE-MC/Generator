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

namespace genie {
namespace utils {
namespace rew   {

  //! Returns a weight to account for a change in hadron mean free path
  double WeightTwkMeanFreePath(double prob_def, double prob_twk, bool interacted);


}  // rew   namespace
}  // utils namespace
}  // genie namespace

#endif // _RW_UTILS_H_
