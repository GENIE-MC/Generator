//____________________________________________________________________________
/*!

\class    genie::FLUKA

\brief    A GENIE interface to the FLUKA hadron transport code.
          Is a concerete implementation of the EventRecordVisitorI interface.
          Note: To use this event generation module you need to obtain FLUKA 
          from its official distribution point and enable it during the GENIE 
          installation.

\ref      More information at: http://www.fluka.org

          FLUKA Authors:
          G.Battistoni, A.Ferrari, P.R.Sala (INFN & Univ. Milano, CERN)

          The FLUKA code is maintained and developed under INFN-CERN agreement
          and copyright 1989-2007. Please cite FLUKA separately if you include 
          this event generation module in your event generation threads.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 13, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FLUKA_H_
#define _FLUKA_H_

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class FLUKA : public EventRecordVisitorI {

public:
  FLUKA();
  FLUKA(string config);
 ~FLUKA();

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
#endif // _FLUKA_H_
