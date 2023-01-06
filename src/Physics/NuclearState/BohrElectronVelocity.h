//____________________________________________________________________________
/*!

\class    genie::BohrElectronVelocity

\brief    It visits the event record & computes a Bohr Velocity for
          initial state electrons bound in coloumb potential.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Brinden Carlson <bcarlson1 \at ufl.edu>
          University of Florida & Fermilab

\created  December 5, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

// #ifndef _BOHR_ELECTRON_VELOCITY_H_
// #define _BOHR_ELECTRON_VELOCITY_H_

#include "Physics/NuclearState/ElectronVelocity.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class BohrElectronVelocity : public ElectronVelocity {
public:
  virtual void InitializeVelocity(Interaction & interaction) const final;
  BohrElectronVelocity();
  BohrElectronVelocity(const string & config);
  ~BohrElectronVelocity();

private:
  float bohr_velocity(int n, int Z) const; //Bohr velocity
  int random_n(int Z) const; //Return random energy level from n_dist
  float random_bohr_velocity(int Z) const; //Generate random n, then calculate velocity from there
};
}