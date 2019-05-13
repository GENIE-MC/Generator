//____________________________________________________________________________
/*!

\class   genie::NtpMCEventRecord

\brief   MINOS-style ntuple record. Each such ntuple record holds a generated
         EventRecord object. Ntuples of this type are intended for feeding
         GENIE events into other applications (for example the GEANT4 based
         MC generation framework of an experiment) if no direct interface
         exists.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 1, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_EVENT_RECORD_H_
#define _NTP_MC_EVENT_RECORD_H_

#include <ostream>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCRecordI.h"

using std::ostream;

namespace genie {

class NtpMCEventRecord;
ostream & operator<< (ostream& stream, const NtpMCEventRecord & rec);

class NtpMCEventRecord : public NtpMCRecordI {

public :
  using NtpMCRecordI::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  NtpMCEventRecord();
  NtpMCEventRecord(const NtpMCEventRecord & ntpmcrec);
  virtual ~NtpMCEventRecord();

  void Fill  (unsigned int ievent, const EventRecord * ev_rec);
  void Copy  (const NtpMCEventRecord & ntpmcrec);
  void Clear (Option_t * opt = "");

  void PrintToStream(ostream & stream) const;
  friend ostream & operator<< (ostream& stream, const NtpMCEventRecord & rec);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  EventRecord * event; ///< event

private:

  void Init (void);

  ClassDef(NtpMCEventRecord, 1)
};

}      // genie namespace

#endif // _NTP_MC_EVENT_RECORD_H_
