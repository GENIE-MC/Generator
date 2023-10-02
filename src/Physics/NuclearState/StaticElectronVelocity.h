//____________________________________________________________________________
/*!

\class    genie::StaticElectronVelocity

\brief    It visits the event record & initializes a static velocity for
          initial state electron.

\author   Brinden Carlson <bcarlson1 \at ufl.edu>
          University of Florida & Fermilab

\created  December 5, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include "Physics/NuclearState/ElectronVelocity.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class StaticElectronVelocity : public ElectronVelocity {
public:
  void InitializeVelocity(Interaction & interaction) const final;

  StaticElectronVelocity();
  StaticElectronVelocity(const string & config);
  ~StaticElectronVelocity();

};
}