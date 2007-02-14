//____________________________________________________________________________
/*!

\class    genie::PrimaryVtxGenerator

\brief    Generates a primary interaction vertex assuming a 'liquid drop' model
          for nuclear targets so to give a 'starting point' for cascading MCs
          (simulating intranuclear effects) that are stepping the interaction
          products out of the nucleus.
          Note that the target is considered to be 'centered' at (0,0,0). When
          running the GENIE's event generation modules using the GENIE MC job
          driver, the driver would shift the vertex at a random point along the
          neutrino direction, in a volume of the specified GEANT/ROOT geometry
          containing the selected target material.
          Is a concrete implementation of the PrimaryVtxGeneratorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 09, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PRIMARY_VERTEX_GENERATOR_H_
#define _PRIMARY_VERTEX_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class PrimaryVtxGenerator : public EventRecordVisitorI {

public :

  PrimaryVtxGenerator();
  PrimaryVtxGenerator(string config);
  ~PrimaryVtxGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _RES_RESONANCE_SELECTOR_H_
