//____________________________________________________________________________
/*!

\class    genie::GlashowResonanceGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_GENERATOR_H_
#define _GLASHOW_RESONANCE_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class GlashowResonanceGenerator : public EventRecordVisitorI {

public :
  GlashowResonanceGenerator();
  GlashowResonanceGenerator(string config);
 ~GlashowResonanceGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

private:

  void PickHitElectron (GHepRecord * event_rec) const;
  void AddResonance    (GHepRecord * event_rec) const;

};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_GENERATOR_H_
