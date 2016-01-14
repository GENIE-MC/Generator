//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::NuXSecFunc

\brief    A set of functors used in gvld_nu_xsec app
          Each functor contains a recipe for extracting a cross-secion out
          of the GENIE simulation inputs.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 17, 2012

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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
// ABC for neutrino cross-section functions used in data/MC comparisons
// performed by code in validation/NuXSec package.
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

  // unique name for each cross-section function
  virtual string Name(void) const 
  { 
    return ""; 
  }
  // initialize
  void Init(GSimFiles * genie_inputs) 
  { 
    fGenieInputs = genie_inputs; 
  }

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
  virtual bool UseTree (EventRecord & /*event*/) { return false; }
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
  bool UseTree (EventRecord & event);
  string Name (void) const;
private:
  int fNuPdg;
  int fTgtPdg;
  int fHitNucPdg;
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
  bool UseTree (EventRecord & event);
  string Name (void) const;
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
  bool UseTree (EventRecord & event);
  string Name (void) const;
private:
  int fNuPdg;
  int fTgtPdg;
  int fPiPdg;
};
//____________________________________________________________________________
// Building cross-section ratio for modes X and Y from a generated event sample 
// treating the generator as a black-box. The cross-section ratio R (E) = 
// sig_{X}(E) / sig_{Y}(E) is calculated from the ratio of generated number
// of events in each energy bin.
// The error envelope is calculated using event-reweighting and considering
// the nuisance params defined in each concrete realization of XSecForModeX.
//
class XSecRatioForModesXY: public NuXSecFunc
{
public:
  XSecRatioForModesXY();
  virtual ~XSecRatioForModesXY();
  virtual TGraphAsymmErrors * ExtractFromEventSample(
     int imodel, double Emin, double Emax, 
     int n, bool inlogE, bool scale_with_E, bool incl_err_band);
//  virtual bool IsCC    (EventRecord & /*event*/) { return false; }
  virtual bool IsModeX (EventRecord & /*event*/) { return false; }
  virtual bool IsModeY (EventRecord & /*event*/) { return false; }
  virtual bool UseTree (EventRecord & /*event*/) { return false; }
protected:
  vector<GSyst_t> fNuisanceParams;    ///<
  GReWeight       fRew;               ///<
//  string          fXSecDirectoryName; ///<
//  string          fModeXSplineName;   ///<
//  string          fModeYSplineName;   ///<
//  double          fXSecScaleFactor;   ///<
};
//____________________________________________________________________________
class CCpi0_CCQE: public XSecRatioForModesXY
{
public:
  CCpi0_CCQE(int nupdg, int tgtpdg);
 ~CCpi0_CCQE();
  bool IsModeX (EventRecord & event);
  bool IsModeY (EventRecord & event);
  bool UseTree (EventRecord & event);
  string Name (void) const;
private:
  int fNuPdg;
  int fTgtPdg;
};
//____________________________________________________________________________
// Inclusive cross-section for scattering off an isoscalar target (n+p)/2
class CCIsoInclXSec: public NuXSecFunc
{
public:
  CCIsoInclXSec(int nupdg);
 ~CCIsoInclXSec();
  TGraphAsymmErrors * ExtractFromEventSample(
     int imodel, double Emin, double Emax, 
     int n, bool inlogE, bool scale_with_E, bool incl_err_band);
  bool UseTree(EventRecord & event);
  string Name (void) const;
private:
  vector<GSyst_t> fNuisanceParams;    ///<
  GReWeight       fRew;               ///<
  int fNuPdg;
};
//____________________________________________________________________________
// Ratio of numubar to numu inclusive cross-section for scattering off an 
// isoscalar target 
class r: public NuXSecFunc
{
public:
  r();
 ~r();
  TGraphAsymmErrors * ExtractFromEventSample(
     int imodel, double Emin, double Emax, 
     int n, bool inlogE, bool scale_with_E, bool incl_err_band);
  bool UseTree(EventRecord & event);
  string Name (void) const;
private:
  vector<GSyst_t> fNuisanceParams;    ///<
  GReWeight       fRew;               ///<
};
//____________________________________________________________________________


} // mc_vs_data namespace
} // genie namepace

#endif  // _NU_XSEC_FUNCTIONS_H_
