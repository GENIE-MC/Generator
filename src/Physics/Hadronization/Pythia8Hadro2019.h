//____________________________________________________________________________
/*!

\class    genie::Pythia8Hadro2019

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
          Queen Mary University of London

\created  October 17, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _PYTHIA8_HADRONIZATION_H_
#define _PYTHIA8_HADRONIZATION_H_

#include "Framework/Conventions/GBuild.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Hadronization/PythiaBaseHadro2019.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Framework/Utils/Pythia8Singleton.h"
#endif

namespace genie {

class GHepParticle;

class Pythia8Hadro2019 : public PythiaBaseHadro2019 {

public:
  Pythia8Hadro2019();
  Pythia8Hadro2019(string config);
  virtual ~Pythia8Hadro2019();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  bool Hadronize (GHepRecord* event) const;

  void CopyOriginalDecayFlags     (void) const;
  void SetDesiredDecayFlags       (void) const;
  void RestoreOriginalDecayFlags  (void) const;

  void LoadConfig (void);
  void Initialize (void);

};

}         // genie namespace

#endif    // _PYTHIA8_HADRONIZATION_H_
