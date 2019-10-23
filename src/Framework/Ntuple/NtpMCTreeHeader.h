//____________________________________________________________________________
/*!

\class   genie::NtpMCTreeHeader

\brief   MINOS-style Ntuple Class to hold an output MC Tree Header

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 1, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_TREE_HEADER_H_
#define _NTP_MC_TREE_HEADER_H_

#include <ostream>

#include <TNamed.h>
#include <TObjString.h>

#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCDTime.h"

using std::string;
using std::ostream;

namespace genie {

class NtpMCTreeHeader;
ostream & operator << (ostream & stream, const NtpMCTreeHeader & hdr);

class NtpMCTreeHeader : public TNamed {

public :
  using TNamed::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  NtpMCTreeHeader();
  NtpMCTreeHeader(const NtpMCTreeHeader & hdr);
  virtual ~NtpMCTreeHeader();

  void Init (void);
  void Copy (const NtpMCTreeHeader & hdr);

  void PrintToStream(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NtpMCTreeHeader & hdr);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breakinsg field data members not prefaced by "f" and mostly lowercase.

  NtpMCFormat_t format;     ///< Event Record format (GENIE support multiple formats)
  TObjString    cvstag;     ///< GENIE CVS Tag (to keep track of GENIE's version)
  NtpMCDTime    datime;     ///< Date and Time that the event ntuple was generated
  Long_t        runnu;      ///< MC Job run number
  TObjString    tune;       ///< GENIE Tune Name
  TObjString    tuneDir;    ///< directory from when tune config came
  TObjString    customDirs; ///< any custom directories

  ClassDef(NtpMCTreeHeader, 3)
};

}      // genie namespace

#endif // _NTP_MC_TREE_HEADER_H_
