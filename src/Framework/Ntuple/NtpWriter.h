//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A utility class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NTP_WRITER_H_
#define _NTP_WRITER_H_

#include <string>

#include "Framework/Ntuple/NtpWriterI.h"
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

class NtpWriter : public NtpWriterI {

public :
  NtpWriter(NtpMCFormat_t fmt = kNFGHEP, Long_t runnu = 0, Long_t runseed = -1);
  virtual ~NtpWriter();

  ///< initialize the ntuple writer
  virtual void Initialize (void) override;

  ///< add event
  virtual void AddEventRecord (int ievent, const EventRecord * ev_rec) override;

  ///< save the event tree
  virtual void Save (void) override;

  ///< get the event tree
  TTree *  EventTree (void) { return fOutTree; }

private:

  virtual void SetDefaultFilename(
    const std::string& filename_prefix = "gntp") override;

  void OpenFile              (string filename);
  void CreateTree            (void);
  void CreateTreeHeader      (void);
  void CreateEventBranch     (void);
  void CreateGHEPEventBranch (void);

  NtpMCFormat_t      fNtpFormat;          ///< enumeration of event formats
  Long_t             fRunNu;              ///< run nu
  Long_t             fRunSeed;            ///< run seed
  TFile *            fOutFile;            ///< output file
  TTree *            fOutTree;            ///< output tree
  TBranch *          fEventBranch;        ///< the generated event branch
  NtpMCEventRecord * fNtpMCEventRecord;   ///<
  NtpMCTreeHeader *  fNtpMCTreeHeader;    ///<
};

}      // genie namespace
#endif // _NTP_WRITER_H_
