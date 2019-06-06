//____________________________________________________________________________
/*!

\class    genie::INukeDeltaPropg

\brief    	 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Oct 01, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INUKE_DELTA_PROPG_H_
#define _INUKE_DELTA_PROPG_H_

#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class INukeDeltaPropg : public EventRecordVisitorI {

public:
  INukeDeltaPropg();
  INukeDeltaPropg(string config);
 ~INukeDeltaPropg();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);

  double fR0;       ///< effective nuclear size param
  double fNR;       ///< param multiplying the nuclear radius, determining how far to track hadrons beyond the "nuclear boundary"
  double fHadStep;  ///< step size for intranuclear hadron transport
};

}      // genie namespace
#endif // _INUKE_DELTA_PROPG_H_
