//____________________________________________________________________________
/*!

\class   genie::exceptions::EVGThreadException

\brief   An exception thrown by EventRecordVisitorI when the normal processing
         sequence has to be disrupted (fast-fwd at the end or step-back)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 27, 2005

*/
//____________________________________________________________________________

#ifndef _EVG_THREAD_EXCEPTION_H_
#define _EVG_THREAD_EXCEPTION_H_

#include <string>
#include <ostream>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class Interaction;

class EVGThreadException {

public :

  EVGThreadException();
  EVGThreadException(const EVGThreadException & exception);
  ~EVGThreadException();

  void   SetReason  (string reason) { fReason   = reason; }
  void   SwitchOnFastForward (void) { fFastFwd  = true;   }
  void   SwitchOnStepBack    (void) { fStepBack = true;   }

  string ShowReason  (void) const { return fReason;   }
  bool   FastForward (void) const { return fFastFwd;  }
  bool   StepBack    (void) const { return fStepBack; }

  void Init  (void);
  void Copy  (const EVGThreadException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
                    ostream & stream, const EVGThreadException & exception);

private:

  bool   fFastFwd;
  bool   fStepBack;
  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _EVG_THREAD_EXCEPTION_H_
