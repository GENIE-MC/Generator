//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include <string>

#include "Framework/Ntuple/NtpWriterI.h"

using namespace genie;

//____________________________________________________________________________
NtpWriterI::~NtpWriterI()
{

}
//____________________________________________________________________________
void NtpWriterI::CustomizeFilename(const std::string& filename)
{
 fOutFilename = filename;
}
//____________________________________________________________________________
void NtpWriterI::CustomizeFilenamePrefix(const std::string& prefix)
{
  this->SetDefaultFilename(prefix);
}
//____________________________________________________________________________
void NtpWriterI::AttachGMCJDriver( const genie::GMCJDriver* /*mc_driver*/ )
{
  // By default, this is a no-op unless a derived class overrides it
}
//____________________________________________________________________________
