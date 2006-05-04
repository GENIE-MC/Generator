//____________________________________________________________________________
/*!

\namespace  genie::utils::interaction

\brief      Interaction Utilities

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

\cpright    Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
            All rights reserved.
            For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INTERACTION_UTILITIES_H_
#define _INTERACTION_UTILITIES_H_

#include "Interaction/Interaction.h"

namespace genie {
namespace utils {

namespace interaction {

  // Get the PDG code of the recoil (final state) nucleon based on the initial
  // state and the process type

  int RecoilNucleonPdgCode    (const Interaction * interaction);
  int QELRecoilNucleonPdgCode (const Interaction * interaction);

  // Provide a more compact way to instantiate interaction objects

  Interaction * GetDis  (int Z, int A, int probe, const TLorentzVector & p4probe, InteractionType_t intype);
  Interaction * GetDis  (int Z, int A, int probe, const InteractionType_t intype);
  Interaction * GetDisCC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetDisCC(int Z, int A, int probe);
  Interaction * GetDisNC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetDisNC(int Z, int A, int probe);

  Interaction * GetQel  (int Z, int A, int probe, const TLorentzVector & p4probe, InteractionType_t intype);
  Interaction * GetQel  (int Z, int A, int probe, const InteractionType_t intype);
  Interaction * GetQelCC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetQelCC(int Z, int A, int probe);
  Interaction * GetQelNC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetQelNC(int Z, int A, int probe);

  Interaction * GetRes  (int Z, int A, int probe, const TLorentzVector & p4probe, InteractionType_t intype);
  Interaction * GetRes  (int Z, int A, int probe, const InteractionType_t intype);
  Interaction * GetResCC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetResCC(int Z, int A, int probe);
  Interaction * GetResNC(int Z, int A, int probe, const TLorentzVector & p4probe);
  Interaction * GetResNC(int Z, int A, int probe);

  Interaction * GetIMD (int Z, int A, int probe, const TLorentzVector & p4probe);

} // interaction namespace

} // utils namespace
} // genie namespace

#endif // _INTERACTION_UTILITIES_H_
