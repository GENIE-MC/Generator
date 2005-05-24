//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2005
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"

namespace genie {

//____________________________________________________________________________
ostream & operator<< (ostream& stream, const XSecSplineList & list)
{
  list.Print(stream);

  return stream;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::fInstance = 0;
//____________________________________________________________________________
XSecSplineList::XSecSplineList()
{
  fInstance    =  0;
  fUseLogE     = true;
  fExtrapolate = false;
  fNKnots      = 100;
  fEmin        =   0.01; // GeV
  fEmax        = 100.00; // GeV
  fEExtrap     = -1.;
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
        const Interaction * interaction, int nknots, double Emin, double Emax)
{
// Build a cross section spline for the input interaction using the input
// cross section algorithm and store in the list.
// For building this specific entry of the spline list, the user is allowed
// to override the list-wide nknots,Emin,Emax
 
  SLOG("XSecSplineList", pINFO) 
             << "Creating cross section spline using the algorithm: " << *alg;

  string key = this->BuildSplineKey(alg,interaction);

  // if any of the nknots,Emin,Emax was not set or its value is not acceptable
  // use the list values
  if (Emin   < 0.) Emin   = this->Emin();
  if (Emax   < 0.) Emin   = this->Emax();
  if (nknots <= 2) nknots = this->NKnots();
  
  assert(Emin < Emax);
  
  double xsec[nknots];
  double E   [nknots];

  bool   first_extarpolation = true;
  double slope = 0, offset = 0;
  
  double dE = 0;  
  if( this->UseLogE() )
       dE = (TMath::Log10(Emax) - TMath::Log10(Emin)) /(nknots-1);
  else dE = (Emax-Emin) /(nknots-1);

  for (int i = 0; i < nknots; i++) {

    if( this->UseLogE() )
          E[i] = TMath::Power(10., TMath::Log10(Emin) + i * dE);
    else  E[i] = Emin + i * dE;

    bool extrapolate = this->Extrapolate() && (E[i] > this->EExtrap()) && (i>1);

    if(extrapolate && first_extarpolation) {
      slope  = (xsec[i-2] - xsec[i-1]) / (E[i-2] - E[i-1]);
      offset = xsec[i-1] - slope * E[i-1];
      first_extarpolation = false;
    }
      
    if(!extrapolate) {
      TLorentzVector p4(0,0,E[i],E[i]);
      interaction->GetInitialStatePtr()->SetProbeP4(p4);

      xsec[i] = alg->XSec(interaction);

      SLOG("XSecSplineList", pINFO)<< "xsec(E = " << E[i] << ") = " << xsec[i];

    } else {
      xsec[i] = offset + slope * E[i];
      
      SLOG("XSecSplineList", pINFO)
            << "xsec(E = " << E[i] << ") = " << xsec[i] << " **extrapolated**";
    }      
  }
  Spline * spline = new Spline(nknots, E, xsec);

  fSplineMap.insert( map<string, Spline *>::value_type(key, spline) );  
}  
//____________________________________________________________________________
void XSecSplineList::SetLogE(bool on)
{
  fUseLogE = on;
}
//____________________________________________________________________________
void XSecSplineList::SetNKnots(int nk)
{
  fNKnots = nk;
  
  if(fNKnots<2) fNKnots = 2; // minimum acceptable number of knots
}
//____________________________________________________________________________
void XSecSplineList::SetMinE(double Ev)
{
  fEmin = Ev;

  if(fEmin<0) fEmin = 0.; 
}
//____________________________________________________________________________
void XSecSplineList::SetMaxE(double Ev)
{
  fEmax = Ev;
  
  if(fEmax<0) fEmax = 0.;
}
//____________________________________________________________________________
void XSecSplineList::SetExtrap(double Ev)
{
  fEExtrap = Ev;

  if(fEExtrap>0) fExtrapolate = true;
  else           fExtrapolate = false;
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
void XSecSplineList::Print(ostream & stream) const
{
  stream << "\n ******************* XSecSplineList *************************";
  stream << "\n [-] Options:";
  stream << "\n  |";
  stream << "\n  |-----o  UseLogE..................." << fUseLogE;
  stream << "\n  |-----o  Extrapolate..............." << fExtrapolate;
  stream << "\n  |-----o  Spline Emin..............." << fNKnots;
  stream << "\n  |-----o  Spline Emax..............." << fEmin;
  stream << "\n  |-----o  Spline NKnots............." << fEmax;
  stream << "\n  |-----o  Extrapolate E............." << fEExtrap;
  stream << "\n  |";  
  stream << "\n [-] Available Splines:";
  stream << "\n  |";

  map<string, Spline *>::const_iterator mapiter;
  for(mapiter = fSplineMap.begin(); mapiter != fSplineMap.end(); ++mapiter) {
    string key = mapiter->first;
    stream << "\n  |-----o  " << key;
  }
  stream << "\n";
}
//___________________________________________________________________________
  
} // genie namespace 
