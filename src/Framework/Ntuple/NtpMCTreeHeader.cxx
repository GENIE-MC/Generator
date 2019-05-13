//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <fstream>
#include <iomanip>

#include <TSystem.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"

using namespace genie;

ClassImp(NtpMCTreeHeader)

using std::ifstream;
using std::endl;
using std::ios;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCTreeHeader & hdr)
  {
     hdr.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCTreeHeader::NtpMCTreeHeader() :
TNamed("header","GENIE output tree header")
{
  this->Init();
}
//____________________________________________________________________________
NtpMCTreeHeader::NtpMCTreeHeader(const NtpMCTreeHeader & hdr) :
TNamed("header","GENIE output tree header")
{
  this->Copy(hdr);
}
//____________________________________________________________________________
NtpMCTreeHeader::~NtpMCTreeHeader()
{

}
//____________________________________________________________________________
void NtpMCTreeHeader::PrintToStream(ostream & stream) const
{
  string sformat = NtpMCFormat::AsString(this->format);
  string scvstag = this->cvstag.GetString().Data();

  stream << "Tree Header Info:"                     << endl
         << "MC run number     -> " << this->runnu  << endl
         << "NtpRecord Format  -> " << sformat      << endl
         << "GENIE CVS Vrs Nu  -> " << scvstag      << endl
         << "File generated at -> " << this->datime << endl;
}
//____________________________________________________________________________
void NtpMCTreeHeader::Copy(const NtpMCTreeHeader & hdr)
{
  this->format = hdr.format;
  this->cvstag.SetString(hdr.cvstag.GetString().Data());
  this->datime.Copy(hdr.datime);
  this->runnu  = hdr.runnu;
}
//____________________________________________________________________________
void NtpMCTreeHeader::Init(void)
{
  string version;
  string genie_path = gSystem->Getenv("GENIE");   
  string filename   = genie_path + "/VERSION";
  bool vrs_file_found = ! (gSystem->AccessPathName(filename.c_str()));
  if (!vrs_file_found) {
   LOG("Ntp", pERROR)
       << "The GENIE version file [" << filename << "] is not accessible";
   version = "NO CVS version number was specified";
  } else {
    ifstream gvinp(filename.c_str(), ios::in);
    gvinp >> version;
    gvinp.close();
  }

  this->format = kNFUndefined;
  this->cvstag.SetString(version.c_str());
  this->datime.Now();
  this->runnu  = 0;
}
//____________________________________________________________________________
