//____________________________________________________________________________
/*!

\class   genie::exceptions::INukeException

\brief   An exception thrown by SimulateHadronState for kinematics problems.
         TwoBodyCollision/Kinematics used a lot, has various failure modes.
         When failure occurs in HAIntranuke, rechoose the fate.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Steve Dytman <dytman \at pitt.edu>
	 Univ. of Pittsburgh         

\created October 10, 2011

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INUKE_EXCEPTION_H_
#define _INUKE_EXCEPTION_H_

#include <string>
#include <ostream>

#include <TMath.h>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class INukeException {

public :
  INukeException();
  INukeException(const INukeException & exception);
 ~INukeException();

  void SetReason(string reason) { fReason = reason; }

  string ShowReason(void) const { return fReason; }

  void Init  (void);
  void Copy  (const INukeException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const INukeException & exception);

private:

  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _INUKE_EXCEPTION_H_
