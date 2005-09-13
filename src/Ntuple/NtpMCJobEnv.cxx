//____________________________________________________________________________
/*!

\class   genie::NtpMCJobEnv

\brief   MINOS-style Ntuple Class to hold a MC Summary Information

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include "Ntuple/NtpMCJobEnv.h"

using std::endl;
using namespace genie;

ClassImp(NtpMCJobEnv)

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCJobEnv & hdr)
  {
     hdr.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCJobEnv::NtpMCJobEnv()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCJobEnv::NtpMCJobEnv(const NtpMCJobEnv & hdr)
{
  this->Copy(hdr);
}
//____________________________________________________________________________
NtpMCJobEnv::~NtpMCJobEnv()
{

}
//____________________________________________________________________________
void NtpMCJobEnv::PrintToStream(ostream & stream) const
{
  stream << "\n*** MC Job Environment:" << endl;
  stream << "$GEVGL    = " << this->gevgl.GetString().Data()    << endl;
  stream << "$GSPLOAD  = " << this->gspload.GetString().Data()  << endl;
  stream << "$GSPSAVE  = " << this->gspsave.GetString().Data()  << endl;
  stream << "$GMSGCONF = " << this->gmsgconf.GetString().Data() << endl;
}
//____________________________________________________________________________
void NtpMCJobEnv::Copy(const NtpMCJobEnv & env)
{
  this->gevgl    = env.gevgl;
  this->gspload  = env.gspload;
  this->gspsave  = env.gspsave;
  this->gmsgconf = env.gmsgconf;
}
//____________________________________________________________________________
void NtpMCJobEnv::Init(void)
{
  this->gevgl.SetString("");
  this->gspload.SetString("");
  this->gspsave.SetString("");
  this->gmsgconf.SetString("");
}
//____________________________________________________________________________

