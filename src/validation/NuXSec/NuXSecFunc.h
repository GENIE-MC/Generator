//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::NuXSecFunc

\brief    A set of functors used in gvld_nu_xsec app
          Each functor contains a recipe for extracting a cross-secion out
          of the GENIE simulation inputs.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Apr 17, 2012

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_XSEC_FUNCTIONS_H_
#define _NU_XSEC_FUNCTIONS_H_

#include <string>
#include <vector>

#include <TGraphAsymmErrors.h>

#include "EVGCore/EventRecord.h"
#include "Utils/GSimFiles.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GSyst.h"

using std::string;
using std::vector;

using namespace genie::rew;

namespace genie {
namespace mc_vs_data {

//____________________________________________________________________________
class NuXSecFunc
{
public:
  virtual ~NuXSecFunc() {}

  virtual TGraphAsymmErrors * ExtractFromEventSample(
     int         /* imodel        */,
     double      /* Emin          */, 
     double      /* Emax          */, 
     int         /* npoints       */, 
     bool        /* inlogE        */, 
     bool        /* scale_with_E  */,
     bool        /* incl_err_band */) 
    { 
      return 0; 
    }

  void Init(GSimFiles * genie_inputs) { fGenieInputs = genie_inputs; }

  static string BuildXSecDirectoryName(int nupdg, int tgtpdg);

protected:

  NuXSecFunc() { fGenieInputs = 0; }
  GSimFiles * fGenieInputs;
};
//____________________________________________________________________________
// Building cross-section for mode X from a generated event sample treating 
// the generator as a black-box. The cross-section sig_{X}(E) is calculated
// as sig_{X}(E) = sig_{CC}(E) * (Nev_{X}(E) / Nev_{CC}(E)).
// The error envelope is calculated using event-reweighting and considering
// the nuisance params defined in each concrete realization of XSecForModeX.
//
class XSecForModeX: public NuXSecFunc
{
public:
  XSecForModeX();
  virtual ~XSecForModeX();
  virtual TGraphAsymmErrors * ExtractFromEventSample(
     int imodel, double Emin, double Emax, 
     int n, bool inlogE, bool scale_with_E, bool incl_err_band);
  virtual bool IsCC    (EventRecord & /*event*/) { return false; }
  virtual bool IsModeX (EventRecord & /*event*/) { return false; }
protected:
  vector<GSyst_t> fNuisanceParams;    ///<
  GReWeight       fRew;               ///<
  string          fXSecDirectoryName; ///<
  string          fModeXSplineName;   ///<
  double          fXSecScaleFactor;   ///<
};
//____________________________________________________________________________
class CCQEXSec: public XSecForModeX
{
public:
  CCQEXSec(int nupdg, int tgtpdg, int hitnucpdg, double xsec_scale=1.);
 ~CCQEXSec();
  bool IsCC    (EventRecord & event);
  bool IsModeX (EventRecord & event);
private:
  int fNuPdg;
  int fTgtPdg;
  int fHitNucPdg;
  int fNTgtNuc;
};
//____________________________________________________________________________
class CCPionXSec: public XSecForModeX
{
public:
  CCPionXSec(
     int nupdg, int tgtpdg, int hitnucpdg, 
     int pip, int npi0, int npim, int np, int nn, double xsec_scale=1.);
 ~CCPionXSec();
  bool IsCC    (EventRecord & event);
  bool IsModeX (EventRecord & event);
private:
  int fNuPdg;
  int fTgtPdg;
  int fHitNucPdg;
  int fNpip;
  int fNpi0; 
  int fNpim; 
  int fNp;
  int fNn;
};
//____________________________________________________________________________
class CohPionXSec: public XSecForModeX
{
public:
  CohPionXSec(int nupdg, int tgtpdg, int pipdg, double xsec_scale=1.);
 ~CohPionXSec();
  bool IsCC    (EventRecord & event);
  bool IsModeX (EventRecord & event);
private:
  int fNuPdg;
  int fTgtPdg;
  int fPiPdg;
};
//____________________________________________________________________________
class CCIsoInclXSec: public NuXSecFunc
{
public:
  CCIsoInclXSec(int nupdg);
 ~CCIsoInclXSec();
  TGraphAsymmErrors * ExtractFromEventSample(
     int imodel, double Emin, double Emax, 
     int n, bool inlogE, bool scale_with_E, bool incl_err_band);
private:
  vector<GSyst_t> fNuisanceParams;    ///<
  GReWeight       fRew;               ///<
  int fNuPdg;
};
//____________________________________________________________________________


} // mc_vs_data namespace
} // genie namepace

#endif  // _NU_XSEC_FUNCTIONS_H_
