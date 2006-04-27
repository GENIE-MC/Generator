//____________________________________________________________________________
/*!

\class   genie::NtpMCRecHeader

\brief   MINOS-style Ntuple Class to hold an MC Event Record Header

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NTP_MC_RECORD_HEADER_H_
#define _NTP_MC_RECORD_HEADER_H_

#include <ostream>

#include <TObject.h>

using std::ostream;

namespace genie {

class NtpMCRecHeader : public TObject {

public :

  NtpMCRecHeader();
  NtpMCRecHeader(const NtpMCRecHeader & hdr);
  virtual ~NtpMCRecHeader();

  void Init (void);
  void Copy (const NtpMCRecHeader & hdr);

  void PrintToStream(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const NtpMCRecHeader & hdr);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  unsigned int  ievent;  ///< Event number

  ClassDef(NtpMCRecHeader, 1)
};

}      // genie namespace

#endif // _NTP_MC_RECORD_HEADER_H_
