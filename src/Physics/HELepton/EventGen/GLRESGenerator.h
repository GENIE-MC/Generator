//____________________________________________________________________________
/*!

\class    genie::GLRESGenerator

\brief    Generator for glashow resonance.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_GENERATOR_H_
#define _GLASHOW_RESONANCE_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

using namespace genie;

namespace genie {

class GLRESGenerator : public EventRecordVisitorI {

public :
  GLRESGenerator();
  GLRESGenerator(string config);
 ~GLRESGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

private:

  void LoadConfig         (void);

  const EventRecordVisitorI * fWDecayer; ///< PYTHIA W decayer

};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_GENERATOR_H_
