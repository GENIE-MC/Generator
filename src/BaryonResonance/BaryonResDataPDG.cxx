//____________________________________________________________________________
/*!

\class    genie::BaryonResDataPDG

\brief    Concrete implementation of the BaryonResDataSetI interface: Its
          configuration registry is loaded from an XML file with PDG baryon
          resonance data and they are served on request.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004
 
*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResDataPDG.h"
#include "BaryonResonance/BaryonResUtils.h"

using std::endl;

using namespace genie;
using namespace genie::utils::res;

//____________________________________________________________________________
BaryonResDataPDG::BaryonResDataPDG() :
BaryonResDataSetI("genie::BaryonResDataPDG")
{
//  fName = "genie::BaryonResDataPDG";
}
//____________________________________________________________________________
BaryonResDataPDG::BaryonResDataPDG(string config) :
BaryonResDataSetI("genie::BaryonResDataPDG", config)
{
//  fName = "genie::BaryonResDataPDG";

//  this->FindConfig();
}
//____________________________________________________________________________
BaryonResDataPDG::~BaryonResDataPDG()
{

}
//____________________________________________________________________________
int BaryonResDataPDG::ResonanceIndex(Resonance_t res) const
{
  int res_idx = 0;
  
  string key = string(AsString(res)) + "-ResIndex";

  if( fConfig->Exists(key) ) fConfig->Get(key, res_idx);

  return res_idx;
}
//____________________________________________________________________________
int BaryonResDataPDG::OrbitalAngularMom(Resonance_t res) const
{
  int L = 0;
  
  string key = string(AsString(res)) + "-L";

  if( fConfig->Exists(key) ) fConfig->Get(key, L);
  
  return L;
}
//____________________________________________________________________________
bool BaryonResDataPDG::IsDeltaResonance(Resonance_t res) const
{
  bool isD = false;
  
  string key = string(AsString(res)) + "-IsDelta";

  if( fConfig->Exists(key) ) fConfig->Get(key, isD);
    
  return isD;
}
//____________________________________________________________________________
bool BaryonResDataPDG::IsNResonance(Resonance_t res) const
{
  bool isN = false;
  
  string key = string(AsString(res)) + "-IsN";

  if( fConfig->Exists(key) ) fConfig->Get(key, isN);
    
  return isN;
}
//____________________________________________________________________________
double BaryonResDataPDG::Mass(Resonance_t res) const
{
  double mass = 0;
  
  string key = string(AsString(res)) + "-Mass";

  if( fConfig->Exists(key) ) fConfig->Get(key, mass);
    
  return mass;
}
//____________________________________________________________________________
double BaryonResDataPDG::Width(Resonance_t res) const
{
  double width = 0;
  
  string key = string(AsString(res)) + "-Width";

  if( fConfig->Exists(key) ) fConfig->Get(key, width);
    
  return width;
}
//____________________________________________________________________________
double BaryonResDataPDG::BreitWignerNorm(Resonance_t res) const
{
  double norm = 0;
  
  string key = string(AsString(res)) + "-BreitWignerNorm";

  if( fConfig->Exists(key) ) fConfig->Get(key, norm);
    
  return norm;
}
//____________________________________________________________________________

