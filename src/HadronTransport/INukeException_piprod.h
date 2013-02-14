//____________________________________________________________________________
/*!

\class   genie::exceptions::INukeException_piprod

\brief   An exception thrown by SimulateHadronState for kinematics problems.
	here, error was in piprod section of Inelastic in HA

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

#ifndef _INUKE_EXCEPTION_H_
#define _INUKE_EXCEPTION_H_

#include <string>
#include <ostream>

#include <TMath.h>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class INukeException_piprod {

public :
  INukeException_piprod();
  INukeException_piprod(const INukeException_piprod & exception);
 ~INukeException_piprod();

  void SetReason(string reason) { fReason = reason; }

  string ShowReason(void) const { return fReason; }

  void Init  (void);
  void Copy  (const INukeException_piprod & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const INukeException_piprod & exception);

private:

  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _INUKE_EXCEPTION_H_
