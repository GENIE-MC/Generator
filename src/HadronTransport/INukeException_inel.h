//____________________________________________________________________________
/*!

\class   genie::exceptions::INukeException

\brief   An exception thrown by SimulateHadronState for kinematics problems.
         TwoBodyCollision/Kinematics used a lot, has various failure modes.
         When failure occurs in HAIntranuke, rechoose the fate.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Steve Dytman <dytman \at pitt.edu>
	 Univ. of Pittsburgh         

\created October 10, 2011

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INUKE_EXCEPTION_INEL_H_
#define _INUKE_EXCEPTION_INEL_H_

#include <string>
#include <ostream>

#include <TMath.h>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class INukeException_inel {

public :
  INukeException_inel();
  INukeException_inel(const INukeException_inel & exception);
 ~INukeException_inel();

  void SetReason(string reason) { fReason = reason; }

  string ShowReason(void) const { return fReason; }

  void Init  (void);
  void Copy  (const INukeException_inel & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const INukeException_inel & exception);

private:

  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _INUKE_EXCEPTION_H_
