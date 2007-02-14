//____________________________________________________________________________
/*!

\class    genie::GIBUU

\brief    A GENIE interface to the GIBUU hadron transport code.
          Is a concerete implementation of the EventRecordVisitorI interface.
          Note: To use this event generation module you need to obtain GIBUU
          from its official distribution point and enable it during the GENIE
          installation.

\ref      http://tp8.physik.uni-giessen.de:8080/GiBUU/

          GIBUU team: 
          U.Model, O.Bub, T.Gaitanos, K.Gallmeister, D.Kalok, M.Kaskulov,
          A.Larionov, T.Leitner, B.Steinmu¼lle, L. Alvarez-Ruso, P. Mu¼hlich 

          Please cite GIBUU separately if you include this event generation
          module in your event generation threads.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  December 13, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GIBUU_H_
#define _GIBUU_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class GIBUU : public EventRecordVisitorI {

public:
  GIBUU();
  GIBUU(string config);
 ~GIBUU();

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
#endif // _GIBUU_H_
