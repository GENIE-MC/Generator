#include <iostream>

#include <TMath.h>

#include "DBUtils/MultiGraph.h"
#include "NuVldGUI/NeuGenFitParams.h"
#include "NuVldGUI/fit_functions.h"
#include "NuVldGUI/NeuGenCards.h"
#include "NuVldGUI/NuVldUserData.h"
#include "Facades/NeuGenWrapper.h"

using std::endl;
using std::cout;

using namespace genie::nuvld;
using namespace genie::nuvld::facades;

//___________________________________________________________________
void update_neugen_config_with_fit_params(double * par)
{
// NeuGEN parameters that can be fitted:
//
// par[0]   --- Axial mass - QEL
// par[1]   --- Axial mass - RES
// par[2]   --- Axial mass - COH
// par[3]   --- QEL - Axial Form Factor FA(Q^2 = 0)
// par[4]   --- QEL - Elastic Scattering Param Eta
// par[5]   --- RES - Rein/Seghal parameter Omega
// par[6]   --- RES - Rein/Seghal parameter Z
// par[7]   --- COH - Nuclear size scale R0
// par[8]   --- COH - Pion scattering amplitude Re/Im
// par[9]   --- KNO Hadronization model, parameter B / v+p
// par[10]  --- KNO Hadronization model, parameter B / v+n
// par[11]  --- KNO Hadronization model, parameter B / vb+p
// par[12]  --- KNO Hadronization model, parameter B / vb+n
// par[13]  --- KNO Hadronization model, parameter A / v+p
// par[14]  --- KNO Hadronization model, parameter A / v+n
// par[15]  --- KNO Hadronization model, parameter A / vb+p
// par[16]  --- KNO Hadronization model, parameter A / vb+n
// par[17]  --- DIS/RES tuning parameter / v+p  / multiplicity = 2
// par[18]  --- DIS/RES tuning parameter / v+p  / multiplicity = 3
// par[19]  --- DIS/RES tuning parameter / v+n  / multiplicity = 2
// par[20]  --- DIS/RES tuning parameter / v+n  / multiplicity = 3
// par[21]  --- DIS/RES tuning parameter / vb+p / multiplicity = 2
// par[22]  --- DIS/RES tuning parameter / vb+p / multiplicity = 3
// par[23]  --- DIS/RES tuning parameter / vb+n / multiplicity = 2
// par[24]  --- DIS/RES tuning parameter / vb+n / multiplicity = 3
// par[25]  --- XSEC NORM
//
// From the NuValidator GUI:
// -- You determine which parameters are actually fitted.
// -- You specify (min, max) limits for the fitted parameters.
// -- You specify the fixed values for the parameters that are not fitted
//
  NeuGenCards * cards = NeuGenCards::Instance();

  cards -> CurrConfig() -> SetMaQel    ( (float) par[0] );
  cards -> CurrConfig() -> SetMaRes    ( (float) par[1] );
  cards -> CurrConfig() -> SetMaCoh    ( (float) par[2] );
  cards -> CurrConfig() -> SetFa0Qel   ( (float) par[3] );
  cards -> CurrConfig() -> SetEtaQel   ( (float) par[4] );
  cards -> CurrConfig() -> SetOmegaRes ( (float) par[5] );
  cards -> CurrConfig() -> SetZRes     ( (float) par[6] );
  cards -> CurrConfig() -> SetR0Coh    ( (float) par[7] );
  cards -> CurrConfig() -> SetREICoh   ( (float) par[8] );
  cards -> CurrConfig() -> SetKnoB     ( e_vp,  (float) par[9] );
  cards -> CurrConfig() -> SetKnoB     ( e_vn,  (float) par[10] );
  cards -> CurrConfig() -> SetKnoB     ( e_vbp, (float) par[11] );
  cards -> CurrConfig() -> SetKnoB     ( e_vbn, (float) par[12] );
  cards -> CurrConfig() -> SetKnoA     ( e_vp,  (float) par[13] );
  cards -> CurrConfig() -> SetKnoA     ( e_vn,  (float) par[14] );
  cards -> CurrConfig() -> SetKnoA     ( e_vbp, (float) par[15] );
  cards -> CurrConfig() -> SetKnoA     ( e_vbn, (float) par[16] );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vp,  (float) par[17] );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vp,  (float) par[18] );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vn,  (float) par[19] );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vn,  (float) par[20] );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vbp, (float) par[21] );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vbp, (float) par[22] );
  cards -> CurrConfig() -> SetDisRes   ( 2, e_vbn, (float) par[23] );
  cards -> CurrConfig() -> SetDisRes   ( 3, e_vbn, (float) par[24] );

//  cout << "CURRENT NeuGEN CONFIG: " << endl;
//  cout << *(cards -> CurrConfig());
}
//___________________________________________________________________
void update_neugen_config(double par, int iparam)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  switch(iparam) {
  case(0) :
      cards ->CurrConfig()->SetMaQel( (float) par );
      break;
  case(1) :
      cards->CurrConfig()->SetMaRes( (float) par );
      break;
  case(2) :
      cards->CurrConfig()->SetMaCoh( (float) par );
      break;
  case(3) :
      cards->CurrConfig()->SetFa0Qel( (float) par );
      break;
  case(4) :
      cards->CurrConfig()->SetEtaQel( (float) par );
      break;
  case(5) :
      cards->CurrConfig()->SetOmegaRes( (float) par );
      break;
  case(6) :
      cards->CurrConfig()->SetZRes( (float) par );
      break;
  case(7) :
      cards->CurrConfig()->SetR0Coh( (float) par );
      break;
  case(8) :
      cards->CurrConfig()->SetREICoh( (float) par );
      break;
  case(9):
      cards->CurrConfig()->SetKnoB( e_vp,  (float) par );
      break;
  case(10):
      cards->CurrConfig()->SetKnoB( e_vn,  (float) par );
      break;
  case(11):
      cards->CurrConfig()->SetKnoB( e_vbp, (float) par );
      break;
  case(12):
      cards->CurrConfig()->SetKnoB( e_vbn, (float) par );
      break;
  case(13):
      cards->CurrConfig()->SetKnoA( e_vp,  (float) par );
      break;
  case(14):
      cards->CurrConfig()->SetKnoA( e_vn,  (float) par );
      break;
  case(15):
      cards->CurrConfig()->SetKnoA( e_vbp, (float) par );
      break;
  case(16):
      cards->CurrConfig()->SetKnoA( e_vbn, (float) par );
      break;
  case(17):
      cards->CurrConfig()->SetDisRes( 2, e_vp,  (float) par );
      break;
  case(18):
      cards->CurrConfig()->SetDisRes( 3, e_vp,  (float) par );
      break;
  case(19):
      cards->CurrConfig()->SetDisRes( 2, e_vn,  (float) par );
      break;
  case(20):
      cards->CurrConfig()->SetDisRes( 3, e_vn,  (float) par );
      break;
  case(21):
      cards->CurrConfig()->SetDisRes( 2, e_vbp, (float) par );
      break;
  case(22):
      cards->CurrConfig()->SetDisRes( 3, e_vbp, (float) par );
      break;
  case(23):
      cards->CurrConfig()->SetDisRes( 2, e_vbn, (float) par );
      break;
  case(24):
      cards->CurrConfig()->SetDisRes( 3, e_vbn, (float) par );
      break;
  default:
      break;
  }

//  cout << "CURRENT NeuGEN CONFIG: " << endl;
//  cout << *(cards -> CurrConfig());
}
//___________________________________________________________________
void fcn_mgr_floating_norm_e(
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag)
{
// Fit chisquare used in MINUIT

  double chisq = 0;
  NuVldUserData * user_data = NuVldUserData::Instance();

  // get selected data as a MultiGraph with errors

  MultiGraph * mgr = user_data->NuXSec()->GetMultiGraph("stat-noE-scale-with-E");

  double x[1];

  // loop over graphs
  int ngr = mgr->NGraphs();
  for(int igr = 0; igr < ngr; igr++) {

     double scale = 1.0;

     // loop over graph points
     int np = mgr->GetGraph(igr)->GetN();
     for (int i=0; i < np; i++) {

        double xp  = mgr->GetGraph(igr)->GetX()[i];
        double yp  = mgr->GetGraph(igr)->GetY()[i];
        double err = mgr->GetGraph(igr)->GetErrorY(i);

        scale = par[kNNGFitParams+igr];

        yp *= scale;
        x[0] = xp;

        double fm = exclusive_xsec_fitfunc_e(x, par);

        double delta = (yp - fm) / err;

        chisq += delta*delta;
     }
     chisq += (1-scale)*(1-scale)/(0.15*0.15);
  }
  f = chisq;
}
//___________________________________________________________________
void fcn_mgr_floating_norm(
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag)
{
// Fit chisquare used in MINUIT

  double chisq = 0;
  NuVldUserData * user_data = NuVldUserData::Instance();

  // get selected data as a MultiGraph with errors

  MultiGraph * mgr = user_data->NuXSec()->GetMultiGraph("stat-noE");

  double x[1];

  // loop over graphs
  int ngr = mgr->NGraphs();
  for(int igr = 0; igr < ngr; igr++) {

     double scale = 1.0;

     // loop over graph points
     int np = mgr->GetGraph(igr)->GetN();
     for (int i=0; i < np; i++) {

        double xp  = mgr->GetGraph(igr)->GetX()[i];
        double yp  = mgr->GetGraph(igr)->GetY()[i];
        double err = mgr->GetGraph(igr)->GetErrorY(i);

        scale = par[kNNGFitParams+igr];

        yp *= scale;
        x[0] = xp;

        double fm = exclusive_xsec_fitfunc(x, par);

        double delta = (yp - fm) / err;

        chisq += delta*delta;
     }
     chisq += (1-scale)*(1-scale)/(0.15*0.15);
  }
  f = chisq;
}
//___________________________________________________________________
double xsec_fitfunc(double * x, double * par)
{
  //-- Update the current NeuGEN config to correspond to fit params

  update_neugen_config_with_fit_params(par);

  //-- Get the input energy

  float energy = (float) x[0];

  //-- Get the interaction & cuts

  NeuGenCards * cards = NeuGenCards::Instance();

  NGInteraction intr = cards->CurrInputs()->GetInteraction();
  NeuGenCuts cuts    = cards->CurrInputs()->GetCuts();

  //-- Reconfigure the neugen wrapper with the updated neugen params

  NeuGenWrapper neugen( cards -> CurrConfig() );

  //-- Ask the neugen wrapper for the x-section

  double xsec = neugen.XSec(energy, &intr, &cuts);

  cout << "XSec (E = " << energy << ") = " << xsec << endl;

  double scale = par[kNNGFitParams];
  xsec = xsec*scale;

  return xsec;
}
//___________________________________________________________________
double xsec_fitfunc_e(double * x, double * par)
{
  float energy = (float) x[0];

  if(energy>0) return xsec_fitfunc(x,par)/energy;
  else         return 0;
}
//___________________________________________________________________
double exclusive_xsec_fitfunc(double * x, double * par)
{
  //-- Update the current NeuGEN config card to correspond to fit params

  update_neugen_config_with_fit_params(par);

  //-- Get the input energy

  float energy = (float) x[0];

  //-- Get the interaction, final state & cuts

  NeuGenCards * cards = NeuGenCards::Instance();

  NGInteraction intr = cards->CurrInputs()->GetInteraction();
  NGFinalState  fs   = cards->CurrInputs()->GetFinalState();
  NeuGenCuts    cuts = cards->CurrInputs()->GetCuts();

  //-- Reconfigure the neugen wrapper with the updated neugen params

  NeuGenWrapper neugen( cards -> CurrConfig() );

  //-- Ask the neugen wrapper for the x-section

  double xsec = neugen.ExclusiveXSec(energy, &intr, &fs, &cuts);

  cout << "XSec (E = " << energy << ") = " << xsec << endl;

  double scale = par[kNNGFitParams];
  xsec = xsec*scale;

  return xsec;
}
//___________________________________________________________________
double exclusive_xsec_fitfunc_e(double * x, double * par)
{
  float energy = (float) x[0];

  if(energy>0) return exclusive_xsec_fitfunc(x,par)/energy;
  else         return 0;
}
//___________________________________________________________________

