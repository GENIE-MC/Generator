//____________________________________________________________________________
/*!

\class   genie::exceptions::CmdLineArgParserException

\brief   An exception thrown by the command line argument parser when command
         line arguments can not be found.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 03, 2005

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
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
