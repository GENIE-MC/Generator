//____________________________________________________________________________
/*!

\namespace genie::intranuke

\brief     INTRANUKE utilities

\author    Jim Dobson <j.dobson07 \at imperial.ac.uk>
           Imperial College London

\created   Mar 03, 2009

\cpright   Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_UTILS_H_
#define _INTRANUKE_UTILS_H_

#include "HadronTransport/INukeHadroFates.h"

class TLorentzVector;

namespace genie {

class GHepRecord;

namespace utils {
namespace intranuke
{
  //! Mean free path
  double MeanFreePath(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double nRpi=0.5, double nRnuc=1.0);

  //! Hadron survival probability
  double ProbSurvival(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double mfp_scale_factor=1.0,
    double nRpi=0.5, double nRnuc=1.0, double NR=3, double R0=1.4);
 
  //! Reconstruct the INTRANUKE/hA model fate for the hadron at position i
  INukeFateHA_t ReconstructHadronFateHA (
    GHepRecord * event, int i, bool hA_mode=false); 

  //! Distance to exit
  double Dist2Exit(
    const TLorentzVector & x4, const TLorentzVector & p4, 
    double A, double NR=3, double R0=1.4);

  //! Distance to exit
  double Dist2ExitMFP(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, 
    double A, double NR=3, double R0=1.4);

}      // intranuke namespace
}      // utils     namespace
}      // genie     namespace

#endif // _INTRANUKE_UTILS_H_
