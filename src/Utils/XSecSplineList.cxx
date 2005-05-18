//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2005
*/
//____________________________________________________________________________

#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"

namespace genie {

//____________________________________________________________________________
XSecSplineList * XSecSplineList::fInstance = 0;
//____________________________________________________________________________
XSecSplineList::XSecSplineList()
{
  fInstance =  0;
}
//____________________________________________________________________________
XSecSplineList::~XSecSplineList()
{
  fInstance = 0;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::Instance()
{
  if(fInstance == 0) {

    static XSecSplineList::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new XSecSplineList;
  }
  return fInstance;
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);

  if(fSplineMap.count(key) != 1) return true;

  return false;
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);

  if ( this->SplineExists(alg,interaction) ) {
  
     map<string, Spline *>::const_iterator iter = fSplineMap.find(key);
     return iter->second;

  } else {     
    SLOG("XSecSplineList", pWARN) << "Couldn't find spline for key = " << key;
    return 0;
  }
  return 0;  
}
//____________________________________________________________________________
void XSecSplineList::CreateSpline(const XSecAlgorithmI * alg,
                   const Interaction * interaction, double Emin, double Emax)
{
  string key = this->BuildSplineKey(alg,interaction);

  double xsec[100];
  double E[100];
  double dE = (Emax-Emin)/99.;

  for (int i = 0; i < 100; i++) {

    E[i] = Emin + i * dE;
    TLorentzVector p4(0,0,E[i],E[i]);

    interaction->GetInitialStatePtr()->SetProbeP4(p4);

    xsec[i] = alg->XSec(interaction);

    SLOG("XSecSplineList", pINFO) << "xsec(E = " << E[i] << ") = " << xsec[i];
  }
  Spline * spline = new Spline(100, E, xsec);

  fSplineMap.insert( map<string, Spline *>::value_type(key, spline) );  
}  
//____________________________________________________________________________
string XSecSplineList::BuildSplineKey(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string alg_name  = alg->Name();
  string param_set = alg->ParamSet();
  string intkey    = interaction->AsString();

  string key = alg_name + "/" + param_set + "/" + intkey;

  return key;
}
//____________________________________________________________________________
  
} // genie namespace 
