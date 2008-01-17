//____________________________________________________________________________
/*!

\class    genie::NEUTCascade

\brief    A GENIE interface to the NEUT hadron transport code developed by
	  the SuperK/K2K experiment (Y.Hayato et al.)

\ref      Y.Hayato, Nucl.Phys.Proc.Suppl.112:171-176, 2002

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 05, 2007

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NEUT_CASCADE_H_
#define _NEUT_CASCADE_H_

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecordVisitorI.h"


#ifdef __GENIE_NEUT_CASCADE_ENABLED__

// C bindings for the NEUT cascade subroutines
//
extern "C" {

}
#endif // __GENIE_GIBUU_ENABLED__


namespace genie {

class NEUTCascade : public EventRecordVisitorI {

public:
  NEUTCascade();
  NEUTCascade(string config);
 ~NEUTCascade();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //-- members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}      // genie namespace
#endif // _NEUT_CASCADE_H_
