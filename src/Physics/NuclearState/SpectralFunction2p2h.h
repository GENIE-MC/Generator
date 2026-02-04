//____________________________________________________________________________
/*!

\class    genie::SpectralFunction2p2h

\brief    Speficif implementation of SecondNucleonEmissionI
          to emit the second nulceon coming from a 2p2h pair
          When GENIE is operating in with EffectiveSF

\author   Afroditi Papadopoulou <apapadop \at mit.edu>
          Massachusetts Institute of Technology - October 04, 2019
          Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  October, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_2P2H_H_
#define _SPECTRAL_FUNCTION_2P2H_H_

#include "Framework/GHEP/GHepParticle.h"

#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/SecondNucleonEmissionI.h"


namespace genie {

class NuclearModelI;

class SpectralFunction2p2h : public SecondNucleonEmissionI {

public :
  SpectralFunction2p2h();
  SpectralFunction2p2h(string config);
 ~SpectralFunction2p2h();
  void LoadConfig (void);

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

};

}      // genie namespace
#endif // _SRC_NUCLEAR_RECOIL_H_
