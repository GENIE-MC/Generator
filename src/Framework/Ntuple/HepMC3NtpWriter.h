//____________________________________________________________________________
/*!

\class    genie::HepMC3NtpWriter

\brief    Writes GENIE events to a HepMC3-format output file

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 7, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HEPMC3_NTP_WRITER_H_
#define _HEPMC3_NTP_WRITER_H_

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_HEPMC3_INTERFACE_ENABLED__

#include <memory>

#include "Framework/Ntuple/NtpWriterI.h"

// Forward-declare needed HepMC3 classes here
namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
  class WriterAscii;
}

namespace genie {

class EventRecord;
class HepMC3Converter;

class HepMC3NtpWriter : public NtpWriterI {

public:

  HepMC3NtpWriter();
  virtual ~HepMC3NtpWriter();

  ///< initialize the ntuple writer
  virtual void Initialize() override;

  ///< add event
  virtual void AddEventRecord( int ievent, const EventRecord* ev_rec ) override;

  ///< save the event tree
  virtual void Save() override;

  virtual void AttachGMCJDriver( const GMCJDriver* mc_driver ) override;

protected:

  virtual void SetDefaultFilename(
    const std::string& filename_prefix = "gntp") override;

  std::shared_ptr< genie::HepMC3Converter > fConverter;
  std::shared_ptr< HepMC3::WriterAscii > fWriter;
};

} // genie namespace

#endif // __GENIE_HEPMC3_INTERFACE_ENABLED__
#endif // _HEPMC3_NTP_WRITER_H_
