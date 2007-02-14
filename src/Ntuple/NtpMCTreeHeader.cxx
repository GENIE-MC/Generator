//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Ntuple/NtpMCTreeHeader.h"

using namespace genie;

ClassImp(NtpMCTreeHeader)

using std::endl;

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
  this->format = kNFUndefined;
  this->cvstag.SetString("NO CVS version number was specified");
  this->datime.Now();
  this->runnu  = 0;
}
//____________________________________________________________________________
