//____________________________________________________________________________
/*!

\class   genie::NtpMCTreeHeader

\brief   MINOS-style Ntuple Class to hold an output MC Tree Header

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

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
  stream << "Tree Header Info:" << endl
         << "NtpRecord Format  -> "
                 << NtpMCFormat::AsString(this->format) << endl
         << "GENIE CVS Vrs Nu  -> ["
                     << this->cvstag.GetString().Data() << "]" << endl
         << "File generated at -> " << this->datime
         << endl;
}
//____________________________________________________________________________
void NtpMCTreeHeader::Copy(const NtpMCTreeHeader & hdr)
{
  this->format = hdr.format;
  this->cvstag.SetString(hdr.cvstag.GetString().Data());
  this->datime.Copy(hdr.datime);
}
//____________________________________________________________________________
void NtpMCTreeHeader::Init(void)
{
  this->format = kNFUndefined;
  this->cvstag.SetString("NO CVS version number was specified");
  this->datime.Now();
}
//____________________________________________________________________________
