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
                        int procmask, bool sumQel, bool sumRes, bool sumDis) :
_cutset   (true), 
_kvid     (kvid), 
_kvmin    (kvmin), 
_kvmax    (kvmax), 
_procmask (procmask),
_sumqel   (sumQel), 
_sumres   (sumRes), 
_sumdis   (sumDis)
{

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
