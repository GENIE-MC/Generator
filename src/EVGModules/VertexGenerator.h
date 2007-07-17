//____________________________________________________________________________
/*!

\class    genie::VertexGenerator

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 16, 2007

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _VERTEX_GENERATOR_H_
#define _VERTEX_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class VertexGenerator : public EventRecordVisitorI {

public :
  VertexGenerator();
  VertexGenerator(string config);
 ~VertexGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void  LoadConfig (void);

  double fR0;
};

}      // genie namespace
#endif // _VERTEX_GENERATOR_H_
