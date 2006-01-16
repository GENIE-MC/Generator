//____________________________________________________________________________
/*!

\class   genie::NtpMCEventRecord

\brief   MINOS-style ntuple record. Each such ntuple record holds a generated
         EventRecord object. Ntuples of this type are intended for feeding
         GENIE events into other applications (for example the GEANT4 based
         MC generation framework of an experiment) if no direct interface
         exists.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _NTP_MC_EVENT_RECORD_H_
#define _NTP_MC_EVENT_RECORD_H_

#include <ostream>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCRecordI.h"

using std::ostream;

namespace genie {

class NtpMCEventRecord : public NtpMCRecordI {

public :

  NtpMCEventRecord();
  NtpMCEventRecord(const NtpMCEventRecord & ntpmcrec);
  virtual ~NtpMCEventRecord();

  void Fill (unsigned int ievent, const EventRecord * ev_rec);
  void Copy (const NtpMCEventRecord & ntpmcrec);

  void PrintToStream(ostream & stream) const;
  friend ostream & operator<< (ostream& stream, const NtpMCEventRecord & rec);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  EventRecord * event; ///< event

private:

  void Init  (void);
  void Clear (void);

  static EventRecord * sevent;

  ClassDef(NtpMCEventRecord, 1)
};

}      // genie namespace

#endif // _NTP_MC_EVENT_RECORD_H_
