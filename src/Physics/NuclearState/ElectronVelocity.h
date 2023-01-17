//____________________________________________________________________________
/*!

\class    genie::ElectronVelocity

\brief  It visits the event record & samples a velocity for
          initial state electrons from a velocity distribution.
        Is a concrete implementation of the EventRecordVisitorI interface.

 \author   Brinden Carlson <bcarlson1 \at ufl.edu>
          University of Florida & Fermilab
  
  \created December 5, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _ELECTRON_VELOCITY_H_
#define _ELECTRON_VELOCITY_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class ElectronVelocity : public EventRecordVisitorI {

public :
  virtual ~ElectronVelocity();
  ElectronVelocity();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const override;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config) override;
  void Configure(string config) override;

  //Make public to use in NuElectronPXsec
  virtual void InitializeVelocity(Interaction & interaction) const = 0; //Give initial velocity

private:
  
protected:
  ElectronVelocity(const string & name, const string & config);
  virtual void LoadConfig (void){;}

};

}      // genie namespace
#endif // _FERMI_MOVER_H_
