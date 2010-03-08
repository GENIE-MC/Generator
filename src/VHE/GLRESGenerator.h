//____________________________________________________________________________
/*!

\class    genie::GLRESGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_GENERATOR_H_
#define _GLASHOW_RESONANCE_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class GLRESGenerator : public EventRecordVisitorI {

public :
  GLRESGenerator();
  GLRESGenerator(string config);
 ~GLRESGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

private:

  void SelectElectronVelocity (GHepRecord * event) const;
  void AddRemnantNucleus      (GHepRecord * event) const;
  void AddResonance           (GHepRecord * event) const;
};

}      // genie namespace

#endif // _GLASHOW_RESONANCE_GENERATOR_H_
