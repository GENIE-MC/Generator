//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A utility class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
class NtpMCEventRecord;
class NtpMCTreeHeader;

class NtpWriter {

public :
  NtpWriter(NtpMCFormat_t fmt = kNFGHEP, Long_t runnu = 0);
 ~NtpWriter();

  void     Initialize       (string filename_prefix="gntp");
  void     AddEventRecord   (int ievent, const EventRecord * ev_rec);
  void     Save             (void);
  TTree *  EventTree        (void) { return fOutTree; }  

private:
  void OpenFile              (string filename_prefix="gntp");
  void CreateTree            (void);
  void CreateTreeHeader      (void);
  void CreateEventBranch     (void);
  void CreateGHEPEventBranch (void);

  NtpMCFormat_t      fNtpFormat;          ///< enumeration of event formats
  Long_t             fRunNu;              ///< run nu
  TFile *            fOutFile;            ///< output file
  TTree *            fOutTree;            ///< output tree
  TBranch *          fEventBranch;        ///< the generated event branch 
  NtpMCEventRecord * fNtpMCEventRecord;   ///< 
  NtpMCTreeHeader *  fNtpMCTreeHeader;    ///<
};

}      // genie namespace
#endif // _NTP_WRITER_H_
