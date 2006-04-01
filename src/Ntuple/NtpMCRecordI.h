//____________________________________________________________________________
/*!

\class   genie::NtpMCRecordI

\brief   MINOS-style base class for ntuple records

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _NTP_MC_RECORD_I_H_
#define _NTP_MC_RECORD_I_H_

#include <TObject.h>

#include "Ntuple/NtpMCRecHeader.h"

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

//  virtual void Init  (void) = 0;
//  virtual void Clear (void) = 0;

ClassDef(NtpMCRecordI, 1)
};

}      // genie namespace

#endif // _NTP_MC_RECORD_I_H_
