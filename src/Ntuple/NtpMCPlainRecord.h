//____________________________________________________________________________
/*!

\class   genie::NtpMCPlainRecord

\brief   MINOS-style ntuple record. Each such ntuple record holds an MC event
         summary and a simplified version of the GHEP event record (as a
         TClones array of NtpMCGHepEntry objects). Ntuples of this type are
         intended for analysis in bare ROOT sessions.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NTP_MC_PLAIN_RECORD_H_
#define _NTP_MC_PLAIN_RECORD_H_

#include <ostream>

#include "Ntuple/NtpMCRecordI.h"
#include "Ntuple/NtpMCSummary.h"

class TClonesArray;

using std::ostream;

namespace genie {

class NtpMCPlainRecord : public NtpMCRecordI {

public :

  NtpMCPlainRecord();
  NtpMCPlainRecord(const NtpMCPlainRecord & ntpmcrec);
  virtual ~NtpMCPlainRecord();

  void Fill  (unsigned int ievent, const EventRecord * ev_rec);
  void Copy  (const NtpMCPlainRecord & ntpmcrec);
  void Clear (Option_t * option = "");

  void PrintToStream(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const NtpMCPlainRecord & rec);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  NtpMCSummary   mc;       ///< MC truth summary
  unsigned int   nentries; ///< number of entries in GHEP array
  TClonesArray * ghep;     ///< simplified GHEP array (NtpMCGHepEntry collection)

private:

  void Init(void);

  static TClonesArray * sghep;

  ClassDef(NtpMCPlainRecord, 1)
};

}      // genie namespace

#endif // _NTP_MC_PLAIN_RECORD_H_
