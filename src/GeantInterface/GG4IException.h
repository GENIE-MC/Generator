//____________________________________________________________________________
/*!

\class   genie::exceptions::GG4IException

\brief   GENIE / GEANT4 Interface Exception

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2005

*/
//____________________________________________________________________________

#ifndef _GG4I_EXCEPTION_H_
#define _GG4I_EXCEPTION_H_

#include <string>
#include <ostream>

using std::string;
using std::ostream;

namespace genie {
namespace exceptions {

class GG4IException {

public :

  GG4IException();
  GG4IException(const GG4IException & exception);
  ~GG4IException();

  void SetReason       (string reason) { fReason       = reason; }
  void SetFileNotFound (bool   flag)   { fFileNotFound = flag;   }
  void SetEmptyTree    (bool   flag)   { fEmptyTree    = flag;   }
  void SetEndOfFile    (bool   flag)   { fEndOfFile    = flag;   }
  void SetFormatError  (bool   flag)   { fFormatError  = flag;   }

  string ShowReason   (void) const { return fReason;       }
  bool   FileNotFound (void) const { return fFileNotFound; }
  bool   EmptyTree    (void) const { return fEmptyTree;    }
  bool   EndOfFile    (void) const { return fEndOfFile;    }
  bool   FormatError  (void) const { return fFormatError;  }

  void Init  (void);
  void Copy  (const GG4IException & exception);
  void Print (ostream & stream) const;

  friend ostream & operator << (
             ostream & stream, const GG4IException & exception);

private:

  bool   fFileNotFound;
  bool   fEmptyTree;
  bool   fEndOfFile;
  bool   fFormatError;
  string fReason;
};

} // exceptions namespace
} // genie namespace

#endif // _GG4I_EXCEPTION_H_
