//____________________________________________________________________________
/*!

\class   genie::NtpWriterI

\brief   A base class for writing GENIE MC Ntuples from GENIE GHEP event
         records

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

\created January 7, 2023

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NTP_WRITER_I_H_
#define _NTP_WRITER_I_H_

#include <string>

namespace genie {

class EventRecord;

/// Defines an interface for GENIE ntuple writers
class NtpWriterI {

public:

  virtual ~NtpWriterI();

  ///< Initialize the ntuple writer at the beginning of a run
  virtual void Initialize(void) = 0;

  ///< Add a new event to the output
  virtual void AddEventRecord(int ievent, const EventRecord* ev_rec) = 0;

  ///< Finalize the output at the end of a run
  virtual void Save(void) = 0;

  ///< Use before Initialize() only if you wish to override the default
  ///< filename, or the default filename prefix
  virtual void CustomizeFilename       (const std::string& filename);
  virtual void CustomizeFilenamePrefix (const std::string& prefix);

protected:

  std::string fOutFilename; ///< output filename

  virtual void SetDefaultFilename(const std::string& filename_prefix) = 0;
};

}      // genie namespace
#endif // _NTP_WRITER_I_H_
