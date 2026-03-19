//____________________________________________________________________________
/*!

\class   genie::exceptions::NievesQELException

\brief   An exception thrown by NievesQELCCPXSec for kinematics problems
         or invalid model configurations (e.g. incompatible nuclear model).
         In case of kinematic failure, xsec is set to 0. Configuration errors
         may prevent further calculation.

\author  Steve Dytman <dytman \at pitt.edu>
	       Univ. of Pittsburgh

         Joe Johnston <jpj13 \at pitt.edu>
	       Univ. of Pittsburgh

\created June 2015

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org         
*/
//____________________________________________________________________________

#ifndef _NIEVES_QEL_EXCEPTION_H_
#define _NIEVES_QEL_EXCEPTION_H_

#include <string>
#include <ostream>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class NievesQELException {

public :
  NievesQELException();
  NievesQELException(const NievesQELException & exception);
 ~NievesQELException();

  void SetReason(const string & reason) { fReason = reason; }

  string ShowReason(void) const { return fReason; }
  void SetIsConfigError(bool v = true) { fIsConfigError = v; }
  bool IsConfigError() const { return fIsConfigError; }

  void Init  (void);
  void Copy  (const NievesQELException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const NievesQELException & exception);

private:

  string fReason;
  bool fIsConfigError = false;
};

}      // exceptions namespace
}      // genie namespace

#endif // _NIEVES_QEL_EXCEPTION_H_
