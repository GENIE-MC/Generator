//____________________________________________________________________________
/*!

\class    genie::Pythia6Hadro2019

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _PYTHIA6_HADRONIZATION_H_
#define _PYTHIA6_HADRONIZATION_H_

#define __GENIE_PYTHIA6_ENABLED__

#include "Framework/Conventions/GBuild.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Hadronization/PythiaBaseHadro2019.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

namespace genie {

class GHepParticle;

class Pythia6Hadro2019 : public PythiaBaseHadro2019 {

public:
  Pythia6Hadro2019();
  Pythia6Hadro2019(string config);
  virtual ~Pythia6Hadro2019();

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

#ifdef __GENIE_PYTHIA6_ENABLED__
  mutable TPythia6 * fPythia;  ///< PYTHIA6 wrapper class
#endif
};

}         // genie namespace

#endif    // _PYTHIA6_HADRONIZATION_H_
