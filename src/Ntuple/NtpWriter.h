//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A utility class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _NTP_WRITER_H_
#define _NTP_WRITER_H_

#include <string>

#include "Ntuple/NtpMCFormat.h"

class TFile;
class TTree;
class TBranch;
class TClonesArray;

using std::string;

namespace genie {

class EventRecord;
class NtpMCPlainRecord;
class NtpMCEventRecord;
class NtpMCTreeHeader;

class NtpWriter {

public :

  NtpWriter(NtpMCFormat_t fmt = kNFEventRecord, Long_t runnu = 0);
  ~NtpWriter();

  void Initialize     (string filename);
  void AddEventRecord (int ievent, const EventRecord * ev_rec);
  void Save           (void);

private:

  void      OpenFile           (string filename);
  void      CreateTree         (void);
  void      CreateTreeHeader   (void);
  TBranch * CreateTreeBranch   (void);
  TBranch * CreatePRTreeBranch (void);
  TBranch * CreateERTreeBranch (void);

  NtpMCFormat_t      fNtpFormat;
  Long_t             fRunNu;
  TFile *            fOutFile;
  TTree *            fOutTree;
  NtpMCPlainRecord * fNtpMCPlainRecord;
  NtpMCEventRecord * fNtpMCEventRecord;
  NtpMCTreeHeader *  fNtpMCTreeHeader;
};

}      // genie namespace

#endif // _NTP_WRITER_H_
