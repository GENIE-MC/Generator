//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenCuts

\brief    NeuGEN's cuts

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include "Facades/NeuGenCuts.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie::nuvld::facades;

ClassImp(NeuGenCuts)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {
    ostream & operator << (ostream & stream, const NeuGenCuts & cuts)
    {
      cuts.Print(stream);
      return stream;
    }
  }  
 }
}    
//____________________________________________________________________________
NeuGenCuts::NeuGenCuts() :
_name("default"),
_cutset(false)
{
}
//____________________________________________________________________________
NeuGenCuts::NeuGenCuts(const char * name) :
_name(name),
_cutset(false)
{

}
//____________________________________________________________________________
NeuGenCuts::NeuGenCuts(NGKineVar_t kvid, float kvmin, float kvmax,
                      bool inclusive, bool incQel, bool incRes, bool incDis) :
_cutset    (true), 
_kvid      (kvid), 
_kvmin     (kvmin), 
_kvmax     (kvmax), 
_inclusive (inclusive), 
_incqel    (incQel), 
_incres    (incRes), 
_incdis    (incDis)
{
  LOG("NuVld",pINFO) << "Updating process mask";

  this->UpdateProcessMask();
}
//____________________________________________________________________________
NeuGenCuts::~NeuGenCuts()
{

}
//____________________________________________________________________________
void NeuGenCuts::SetCut(NGKineVar_t kvid, float kvmin, float kvmax)
{
  _cutset = true;
  _kvid   = kvid;
  _kvmin  = kvmin;
  _kvmax  = kvmax;
}
//____________________________________________________________________________
void NeuGenCuts::SetInclusive(bool on) 
{ 
  _inclusive = on; 
}
//____________________________________________________________________________
void NeuGenCuts::SetIncludeQel(bool on) 
{ 
  _incqel = on; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::SetIncludeRes(bool on) 
{ 
  _incres = on; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::SetIncludeDis(bool on) 
{ 
  _incdis = on; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::UpdateProcessMask(void)
{
  int qel, dis, res;

  (_incqel) ? qel = 1 : qel = 0;
  (_incres) ? res = 1 : res = 0;
  (_incdis) ? dis = 1 : dis = 0;

  _procmask = 15 - 4*dis - 2*res - qel;
}
//____________________________________________________________________________
bool NeuGenCuts::SumQel(void) const 
{ 
  return (_incqel && _inclusive);   
}
//____________________________________________________________________________
bool NeuGenCuts::SumRes(void) const 
{ 
  return (_incres && _inclusive);   
}
//____________________________________________________________________________
bool NeuGenCuts::SumDis(void) const 
{ 
  return (_incdis && _inclusive); 
}
//____________________________________________________________________________
void NeuGenCuts::Print(ostream & stream) const
{
  stream << endl;
  stream << "cut-set:.........." << _cutset                     << endl;
  stream << "Kvid:............." << NGKineVar::AsString(_kvid)  << endl;
  stream << "min:.............." << this->KVMin()               << endl;
  stream << "max:.............." << this->KVMax()               << endl;
  stream << "procmask:........." << this->ProcMask()            << endl;
  stream << "inclusive:........" << this->Inclusive()           << endl;
  stream << "include(qel):....." << this->IncludeQel()          << endl;
  stream << "include(res):....." << this->IncludeRes()          << endl;
  stream << "include(dis):....." << this->IncludeDis()          << endl;
  stream << "sum(qel):........." << this->SumQel()              << endl;
  stream << "sum(res):........." << this->SumRes()              << endl;
  stream << "sum(dis):........." << this->SumDis()              << endl;
}
//____________________________________________________________________________
