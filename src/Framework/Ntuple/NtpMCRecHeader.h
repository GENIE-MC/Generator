//____________________________________________________________________________
/*!

\class   genie::NtpMCRecHeader

\brief   MINOS-style Ntuple Class to hold an MC Event Record Header

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 1, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_RECORD_HEADER_H_
#define _NTP_MC_RECORD_HEADER_H_

#include <ostream>

#include <TObject.h>

using std::ostream;

namespace genie {

class NtpMCRecHeader;
ostream & operator << (ostream & stream, const NtpMCRecHeader & hdr);

class NtpMCRecHeader : public TObject {

public :
  using TObject::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

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
