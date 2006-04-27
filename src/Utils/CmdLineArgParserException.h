//____________________________________________________________________________
/*!

\class   genie::exceptions::CmdLineArgParserException

\brief   An exception thrown by the command line argument parser when command
         line arguments can not be found.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2005

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _CLAP_EXCEPTION_H_
#define _CLAP_EXCEPTION_H_

#include <string>
#include <ostream>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class CmdLineArgParserException {

public :

  CmdLineArgParserException();
  CmdLineArgParserException(const CmdLineArgParserException & exception);
  ~CmdLineArgParserException();

  void SetReason  (string reason) { fReason   = reason; }
  void SwitchOnArgNotFound (void) { fArgFound = false;  }

  string ShowReason    (void) const { return fReason;   }
  bool   ArgumentFound (void) const { return fArgFound; }

  void Init  (void);
  void Copy  (const CmdLineArgParserException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
       ostream & stream, const CmdLineArgParserException & exception);

private:

  bool   fArgFound;
  string fReason;
};

} // exceptions namespace
} // genie namespace

#endif // _CLAP_EXCEPTION_H_
