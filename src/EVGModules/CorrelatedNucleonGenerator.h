//____________________________________________________________________________
/*!

\class    genie::CorrelatedNucleonGenerator

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 16, 2007

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CORRELATED_NUCLEON_GENERATOR_H_
#define _CORELLATED_NUCLEON_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class CorrelatedNucleonGenerator : public EventRecordVisitorI {

public :
  CorrelatedNucleonGenerator();
  CorrelatedNucleonGenerator(string config);
 ~CorrelatedNucleonGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void   LoadConfig (void);

  bool   fSimulateCorrelN;
  double fMomentumThr;
};

}      // genie namespace

#endif // _CORRELATED_NUCLEON_GENERATOR_H_
