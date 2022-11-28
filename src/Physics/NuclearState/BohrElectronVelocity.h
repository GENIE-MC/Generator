//____________________________________________________________________________
/*!

\class    genie::ElectronVelocity

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  October 08, 2004

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

  //void Configure(const Registry & config);
  //void Configure(string config);

  //BohrElectronVelocity();
  BohrElectronVelocity(const string & config);
  ~BohrElectronVelocity();

private:
  float bohr_velocity(int n, int Z) const; //Bohr velocity
  float random_n(int Z) const; //Return random energy level from n_dist
  float random_bohr_velocity(int Z) const; //Generate random n, then calculate velocity from there

  //const int static fMaxElectrons; //Max number of electrons 
  //static std::array<int,6> fnprobs; //Probability dist.

};
}