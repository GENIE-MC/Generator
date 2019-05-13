//____________________________________________________________________________
/*!

\class    genie::VertexGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 16, 2007

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _VERTEX_GENERATOR_H_
#define _VERTEX_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"

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

  //-- Generate the vertex position
  //   public so other classes can reuse this code to generate a position
  TVector3 GenerateVertex(const Interaction * in,double A) const;

private:
  void  LoadConfig (void);

  int    fVtxGenMethod; ///< vtx generation method (0: uniform, 1: according to nuclear density [def])
  double fR0;           ///< parameter controlling nuclear sizes
};

}      // genie namespace
#endif // _VERTEX_GENERATOR_H_
