//____________________________________________________________________________
/*!

\class   genie::NtpMCTreeHeader

\brief   MINOS-style Ntuple Class to hold an output MC Tree Header

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _NTP_MC_TREE_HEADER_H_
#define _NTP_MC_TREE_HEADER_H_

#include <ostream>

#include <TNamed.h>

#include "Ntuple/NtpMCFormat.h"

using std::ostream;

namespace genie {

class NtpMCTreeHeader : public TNamed {

public :

  NtpMCTreeHeader();
  NtpMCTreeHeader(const NtpMCTreeHeader & hdr);
  virtual ~NtpMCTreeHeader();

  void Init (void);
  void Copy (const NtpMCTreeHeader & hdr);

  void PrintToStream(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NtpMCTreeHeader & hdr);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breakinsg field data members not prefaced by "f" and mostly lowercase.

  NtpMCFormat_t format;  ///< Event Record format (GENIE support multiple formats)

  ClassDef(NtpMCTreeHeader, 1)
};

}      // genie namespace

#endif // _NTP_MC_TREE_HEADER_H_
