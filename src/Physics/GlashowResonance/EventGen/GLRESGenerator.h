//____________________________________________________________________________
/*!

\class    genie::GLRESGenerator

\brief    Glashow resonance event generator

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_GENERATOR_H_
#define _GLASHOW_RESONANCE_GENERATOR_H_

#include <TPythia6.h>

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class GLRESGenerator : public EventRecordVisitorI {

public :
  GLRESGenerator();
  GLRESGenerator(string config);
 ~GLRESGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig(void);

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_GENERATOR_H_
