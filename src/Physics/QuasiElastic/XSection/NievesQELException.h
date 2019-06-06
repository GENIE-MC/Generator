//____________________________________________________________________________
/*!

\class   genie::exceptions::NievesQELException

\brief   An exception thrown by NievesQELCCPXSec for kinematics problems.
         When failure occurs, set xsec = 0.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Steve Dytman <dytman \at pitt.edu>
	 Univ. of Pittsburgh         
	 
	 Joe Johnston <jpj13 \at pitt.edu>
	 Univ. of Pittsburgh

\created June 2015

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NIEVES_QEL_EXCEPTION_H_
#define _NIEVES_QEL_EXCEPTION_H_

#include <string>
#include <ostream>

#include <TMath.h>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class NievesQELException {

public :
  NievesQELException();
  NievesQELException(const NievesQELException & exception);
 ~NievesQELException();

  void SetReason(string reason) { fReason = reason; }

  string ShowReason(void) const { return fReason; }

  void Init  (void);
  void Copy  (const NievesQELException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const NievesQELException & exception);

private:

  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _NIEVES_QEL_EXCEPTION_H_
