//____________________________________________________________________________
/*!

\class    genie::CacheBranchNtp

\brief    A simple cache branch storing the cached data in a TNtupleD

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 26, 2004

*/
//____________________________________________________________________________

#include <TNtupleD.h>

#include "Utils/CacheBranchNtp.h"

using namespace genie;

ClassImp(CacheBranchNtp);

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const CacheBranchNtp & cbntp)
  {
     cbntp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
CacheBranchNtp::CacheBranchNtp(void) :
CacheBranchI()
{
  this->Init();
}
//____________________________________________________________________________
CacheBranchNtp::CacheBranchNtp(string name, string branch_def) :
CacheBranchI()
{
  this->Init();
  this->CreateNtuple(name, branch_def);
}
//____________________________________________________________________________
CacheBranchNtp::~CacheBranchNtp()
{
  this->CleanUp();
}
//____________________________________________________________________________
void CacheBranchNtp::Init(void)
{
  fNtp = 0;
}
//____________________________________________________________________________
void CacheBranchNtp::CleanUp(void)
{
  if(fNtp) delete fNtp;
}
//____________________________________________________________________________
void CacheBranchNtp::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//____________________________________________________________________________
void CacheBranchNtp::CreateNtuple(string name, string branch_def)
{
  this->Reset();

  fNtp = new TNtupleD(name.c_str(), "CacheBranchNtp", branch_def.c_str());
  fNtp->SetDirectory(0);
  fNtp->SetCircular(1600000);
}
//____________________________________________________________________________
void CacheBranchNtp::Print(ostream & stream) const
{
  if(fNtp) {
    stream << "type: [CacheBranchNtp]  - nentries: " << fNtp->GetEntries();
  } else {
    stream << " *** NULL ***";
  }
}
//____________________________________________________________________________
TNtupleD * CacheBranchNtp::operator () (void) const
{
  return this->Ntuple();
}
//____________________________________________________________________________
