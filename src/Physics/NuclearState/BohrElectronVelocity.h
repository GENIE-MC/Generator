//____________________________________________________________________________
/*!

\class    genie::BohrElectronVelocity

\brief    It visits the event record & computes a Bohr Velocity for
          initial state electrons bound in coloumb potential.
          Is a concrete implementation of the ElectronVelocity interface.

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
  void InitializeVelocity(Interaction & interaction) const final;
  BohrElectronVelocity();
  BohrElectronVelocity(const string & config);
  ~BohrElectronVelocity();

private:
  double reduced_mass(double nucleus_mass) const; //Reduced mass calculation
  double bohr_velocity(unsigned int fn, unsigned int fZ, double mu) const; //Bohr velocity
  unsigned int random_n(unsigned int fZ) const; //Return random energy level from n_dist
  double random_bohr_velocity(unsigned int fZ, double mu) const; //Generate random n, then calculate velocity from there
  TVector3 randomize_direction_sphere(double v) const;
};
}