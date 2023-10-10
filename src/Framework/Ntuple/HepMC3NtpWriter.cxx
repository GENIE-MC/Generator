#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_HEPMC3_INTERFACE_ENABLED__
//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE includes
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/Ntuple/HepMC3NtpWriter.h"
#include "Framework/Utils/RunOpt.h"

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"

//____________________________________________________________________________
genie::HepMC3NtpWriter::HepMC3NtpWriter()
{
}
//____________________________________________________________________________
genie::HepMC3NtpWriter::~HepMC3NtpWriter()
{

}
//____________________________________________________________________________
void genie::HepMC3NtpWriter::SetDefaultFilename(
  const std::string& filename_prefix )
{
  fOutFilename = filename_prefix + ".hepmc3";
}
//____________________________________________________________________________
void genie::HepMC3NtpWriter::Initialize()
{
  fWriter = std::make_shared< HepMC3::WriterAscii >( fOutFilename );
  fConverter = std::make_shared< genie::HepMC3Converter >();
}
//____________________________________________________________________________
void genie::HepMC3NtpWriter::Save()
{

}
//____________________________________________________________________________
void genie::HepMC3NtpWriter::AddEventRecord( int ievent,
  const EventRecord* ev_rec )
{
  auto evt = fConverter->ConvertToHepMC3( *ev_rec );
  evt->set_event_number( ievent );
  fWriter->write_event( *evt );
}
#endif //__GENIE_HEPMC3_INTERFACE_ENABLED__
