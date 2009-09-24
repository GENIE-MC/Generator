
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

#include <TLorentzVector.h>

#include "EVGCore/EventRecord.h"
#include "ReWeight/GSyst.h"

namespace genie {
namespace utils {
namespace rew   {

  // Returns a weight to account for a change in hadron mean free path
  double MeanFreePathWeight(
    int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
    double mfp_scale_factor, bool interacted,
    double nRpi=0.5, double nRnuc=1.0, double NR=3, double R0=1.4);
  double MeanFreePathWeight(
      double prob_def, double prob_twk, bool interacted);

  // Return the fraction of the hadron rescatering fate described by the input
  // systematic enumeration at the input hadron kinetic energy
  double FateFraction(genie::rew::GSyst_t syst, double kinE, double frac_scale_factor=1.);

  // Check whether the input event is hadronized by AGKY
  bool  HadronizedByAGKY(const EventRecord & event);

  // Check whether the input event is hadronized by AGKY/PYTHIA
  bool  HadronizedByAGKYPythia(const EventRecord & event);
  
  // Compute the hadronic system 4-momentum @ LAB
  TLorentzVector Hadronic4pLAB(const EventRecord & event);

  //
  double AGKYWeight(int pdgc, double xF, double pT2);

}  // rew   namespace
}  // utils namespace
}  // genie namespace

#endif // _RW_UTILS_H_
