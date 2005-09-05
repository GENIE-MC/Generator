//____________________________________________________________________________
/*!

\class   genie::NtpMCRecHeader

\brief   MINOS-style Ntuple Class to hold a MC Summary Information

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

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
NtpMCRecHeader::NtpMCRecHeader()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCRecHeader::NtpMCRecHeader(const NtpMCRecHeader & hdr)
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

