//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Ntuple/NtpMCRecHeader.h"

using namespace genie;

ClassImp(NtpMCRecHeader)

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCRecHeader & hdr)
  {
     hdr.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCRecHeader::NtpMCRecHeader() :
TObject()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCRecHeader::NtpMCRecHeader(const NtpMCRecHeader & hdr) :
TObject()
{
  this->Copy(hdr);
}
//____________________________________________________________________________
NtpMCRecHeader::~NtpMCRecHeader()
{

}
//____________________________________________________________________________
void NtpMCRecHeader::PrintToStream(ostream & stream) const
{
  stream << "\n\n*** Event #: " << this->ievent;
}
//____________________________________________________________________________
void NtpMCRecHeader::Copy(const NtpMCRecHeader & hdr)
{
  this->ievent = hdr.ievent;
}
//____________________________________________________________________________
void NtpMCRecHeader::Init(void)
{
  this->ievent = 0;
}
//____________________________________________________________________________

