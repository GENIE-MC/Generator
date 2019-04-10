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

#include "Framework/Ntuple/NtpMCRecHeader.h"

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

