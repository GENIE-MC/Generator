//____________________________________________________________________________
/*!

\class   genie::NtpMCRecordI

\brief   MINOS-style base class for ntuple records

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 1, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_RECORD_I_H_
#define _NTP_MC_RECORD_I_H_

#include <TObject.h>

#include "Framework/Ntuple/NtpMCRecHeader.h"

namespace genie {

class EventRecord;

class NtpMCRecordI : public TObject {

public :
  virtual ~NtpMCRecordI();

  virtual void Fill(unsigned int ievent, const EventRecord * ev_rec) = 0;

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  NtpMCRecHeader hdr;   ///< record header

protected:
  NtpMCRecordI();

ClassDef(NtpMCRecordI, 1)
};

}      // genie namespace
#endif // _NTP_MC_RECORD_I_H_
