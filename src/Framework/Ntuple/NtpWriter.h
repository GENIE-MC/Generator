//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A utility class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created October 1, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NTP_WRITER_H_
#define _NTP_WRITER_H_

#include <string>

#include "Framework/Ntuple/NtpMCFormat.h"

class TFile;
class TTree;
class TBranch;
class TClonesArray;

using std::string;

namespace genie {

class EventRecord;
class NtpMCEventRecord;
class NtpMCTreeHeader;

class NtpWriter {

public :
  NtpWriter(NtpMCFormat_t fmt = kNFGHEP, Long_t runnu = 0, Long_t runseed = -1);
 ~NtpWriter();

  ///< initialize the ntuple writer
  void Initialize (void);

  ///< add event
  void AddEventRecord (int ievent, const EventRecord * ev_rec);

  ///< save the event tree
  void Save (void);

  ///< get the even tree
  TTree *  EventTree (void) { return fOutTree; }

  ///< use before Initialize() only if you wish to override the default
  ///< filename, or the default filename prefix
  void CustomizeFilename       (string filename);
  void CustomizeFilenamePrefix (string prefix);

private:

  void SetDefaultFilename    (string filename_prefix="gntp");
  void OpenFile              (string filename);
  void CreateTree            (void);
  void CreateTreeHeader      (void);
  void CreateEventBranch     (void);
  void CreateGHEPEventBranch (void);

  NtpMCFormat_t      fNtpFormat;          ///< enumeration of event formats
  Long_t             fRunNu;              ///< run nu
  Long_t             fRunSeed;            ///< run seed
  string             fOutFilename;        ///< output filename
  TFile *            fOutFile;            ///< output file
  TTree *            fOutTree;            ///< output tree
  TBranch *          fEventBranch;        ///< the generated event branch
  NtpMCEventRecord * fNtpMCEventRecord;   ///<
  NtpMCTreeHeader *  fNtpMCTreeHeader;    ///<
};

}      // genie namespace
#endif // _NTP_WRITER_H_
