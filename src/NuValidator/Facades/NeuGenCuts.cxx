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
                                      bool sumQel, bool sumRes, bool sumDis) :
_cutset   (true), 
_kvid     (kvid), 
_kvmin    (kvmin), 
_kvmax    (kvmax), 
_sumqel   (sumQel), 
_sumres   (sumRes), 
_sumdis   (sumDis)
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
void NeuGenCuts::SetSumQel(bool sum) 
{ 
  _sumqel = sum; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::SetSumRes(bool sum) 
{ 
  _sumres = sum; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::SetSumDis(bool sum) 
{ 
  _sumdis = sum; 
  this->UpdateProcessMask();
}
//____________________________________________________________________________
void NeuGenCuts::UpdateProcessMask(void)
{
  int qel, dis, res;

  (_sumqel) ? qel = 0 : qel = 1;
  (_sumres) ? res = 0 : res = 1;
  (_sumdis) ? dis = 0 : dis = 1;

  _procmask = qel + 2 * res + 4 * dis;
}
//____________________________________________________________________________
void NeuGenCuts::Print(ostream & stream) const
{
  stream << "CutSet:..........." << _cutset                     << endl;
  stream << "Kvid:............." << NGKineVar::AsString(_kvid)  << endl;
  stream << "min:.............." << _kvmin                      << endl;
  stream << "max:.............." << _kvmax                      << endl;
  stream << "procmask:........." << _procmask                   << endl;
  stream << "sum(qel):........." << _sumqel                     << endl;
  stream << "sum(res):........." << _sumres                     << endl;
  stream << "sum(dis):........." << _sumdis                     << endl;
}
//____________________________________________________________________________
