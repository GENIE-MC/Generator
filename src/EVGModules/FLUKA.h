//____________________________________________________________________________
/*!

\class   genie::FLUKA

\brief   A GENIE interface to the FLUKA hadron transport code.

         Is a concerete implementation of the EventRecordVisitorI interface.

         Note: the FLUKA code is *not included* in your GENIE installation.
         You need to obtain FLUKA from its official distribution point.

\ref     More information at: http://www.fluka.org

         FLUKA Authors:
         G.Battistoni, A.Ferrari, P.R.Sala (INFN & Univ. Milano, CERN)

         The FLUKA code is maintained and developed under INFN-CERN agreement
         and copyright 1989-2005.

         Please cite FLUKA separately if you include this event generation
         module in your event generation threads.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 13, 2005

*/
//____________________________________________________________________________

#ifndef _FLUKA_H_
#define _FLUKA_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class FLUKA : public EventRecordVisitorI {

public :

  FLUKA();
  FLUKA(string config);
  ~FLUKA();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _FLUKA_H_
