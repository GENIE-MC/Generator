//____________________________________________________________________________
/*!

\class    genie::GLRESP6Generator

\brief    Glashow resonance event generator

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_P6_GENERATOR_H_
#define _GLASHOW_RESONANCE_P6_GENERATOR_H_

#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class GLRESP6Generator : public EventRecordVisitorI {

public :
  GLRESP6Generator();
  GLRESP6Generator(string config);
 ~GLRESP6Generator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig(void);

#ifdef __GENIE_PYTHIA6_ENABLED__
  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
#endif
};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_P6_GENERATOR_H_
