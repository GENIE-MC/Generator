//____________________________________________________________________________
/*!

\class   genie::exceptions::EVGThreadException

\brief   An exception thrown by EventRecordVisitorI when the normal processing
         sequence has to be disrupted (fast-fwd at the end or step-back)

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 27, 2005

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _EVG_THREAD_EXCEPTION_H_
#define _EVG_THREAD_EXCEPTION_H_

#include <string>
#include <ostream>

#include <TMath.h>

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

  void   SetReason  (string reason) { fReason     = reason;          }
  void   SwitchOnFastForward (void) { fFastFwd    = true;            }
  void   SwitchOnStepBack    (void) { fStepBack   = true;            }
  void   SetReturnStep (int s)      { fReturnStep = TMath::Max(0,s); }

  string ShowReason  (void) const { return fReason;     }
  bool   FastForward (void) const { return fFastFwd;    }
  bool   StepBack    (void) const { return fStepBack;   }
  int    ReturnStep  (void) const { return fReturnStep; }

  void Init  (void);
  void Copy  (const EVGThreadException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const EVGThreadException & exception);

private:

  bool   fFastFwd;
  bool   fStepBack;
  int    fReturnStep;
  string fReason;
};

}      // exceptions namespace
}      // genie namespace

#endif // _EVG_THREAD_EXCEPTION_H_
