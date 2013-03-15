//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

/*
#include <cassert>

#include <TGTextEdit.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TMath.h>

#include "Facades/NeuGenWrapper.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include  "ValidationTools/NuVld/GuiFitKernel.h"
#include  "ValidationTools/NuVld/NuVldUserData.h"
#include  "ValidationTools/NuVld/NeuGenCards.h"
#include  "ValidationTools/NuVld/NeuGenFitParams.h"
#include  "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include  "ValidationTools/NuVld/fit_functions.h"
#include "Utils/StringUtils.h"

using genie::RandomGen;
using genie::Spline;
using namespace genie::nuvld;
using namespace genie::nuvld::facades;
using namespace genie::utils::str;

//______________________________________________________________________________
GuiFitKernel::GuiFitKernel()
{
  fFunc1d = 0;

  this->Reset();
}
//______________________________________________________________________________
GuiFitKernel::~GuiFitKernel()
{

}
//______________________________________________________________________________
void GuiFitKernel::Reset(void)
{
  fScaleWithE = false;

  fXmin = 0.;
  fXmax = 0.;

  if( fFunc1d ) delete fFunc1d;

  fFunc1d = 0;
  fNGFP   = 0;

  //tmp names
  lowb    = 0;
  highb   = 0;
  chisq1d = 0;
  chisq2d = 0;
}
//______________________________________________________________________________
void GuiFitKernel::SetFitRange(float xmin, float xmax)
{
  fXmin = xmin;
  fXmax = xmax;
}
//______________________________________________________________________________
void GuiFitKernel::PrintConfig(void)
{
  int np = fNGFP->NFittedParams();

  //-- send to msg service stream

  LOG("NuVld", pINFO) << "N not-fixed fit parameters: " << np;
  LOG("NuVld", pINFO) << "Energy range = [" << fXmin << ", " << fXmax << "]";
  LOG("NuVld", pINFO) << "Scale with energy = [" << fScaleWithE;

  //-- send to GUI's session log TGTextEdit

  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog->Log()->AddLine( Concat("N not-fixed fit parameters: ", np) );
  syslog->Log()->AddLine( Concat("(free var) E-min: ", fXmin)        );
  syslog->Log()->AddLine( Concat("(free var) E-max: ", fXmax)        );
  syslog->Log()->AddLine( Concat("scale with energy: ", fScaleWithE) );
}
//______________________________________________________________________________
void GuiFitKernel::DoSimpleFit(bool fit_norm)
{
  GuiSysLogSingleton * syslog  = GuiSysLogSingleton::Instance();
  NuVldUserData * user_data = NuVldUserData::Instance();

  // free parameter range
  LOG("NuVld", pDEBUG) << "Energy = [" << fXmin << ", " << fXmax << "]";

  // find out if we are fitting xsecs for an exclusive or inclusive channel
  NeuGenCards * cards = NeuGenCards::Instance();
  bool is_inclusive  = cards->CurrInputs()->Inclusive();

  // get a fit function (TF1 object)

  if( fFunc1d ) delete fFunc1d;

  fFunc1d = 0;

  if( fScaleWithE ) {

     if(is_inclusive) {
        LOG("NuVld", pDEBUG) << "Using fit function: xsec_fitfunc_e";
        syslog->Log()->AddLine( "Using fit function: xsec_fitfunc_e" );

        fFunc1d = new TF1("fitfunc", xsec_fitfunc_e, fXmin, fXmax, 1+kNNGFitParams);

     } else {
        LOG("NuVld", pDEBUG) << "Using fit function: exclusive_xsec_fitfunc_e";
        syslog->Log()->AddLine( "Using fit function: exclusive_xsec_fitfunc_e" );

        fFunc1d = new TF1("fitfunc", exclusive_xsec_fitfunc_e, fXmin, fXmax, 1+kNNGFitParams);
     }

  } else {
     if(is_inclusive) {
        LOG("NuVld", pDEBUG) << "Using fit function: xsec_fitfunc";
        syslog->Log()->AddLine( "Using fit function: xsec_fitfunc" );

        fFunc1d = new TF1("fitfunc", xsec_fitfunc, fXmin, fXmax, 1+kNNGFitParams);

     } else {
        LOG("NuVld", pDEBUG) << "Using fit function: exclusive_xsec_fitfunc";
        syslog->Log()->AddLine( "Using fit function: exclusive_xsec_fitfunc" );

        fFunc1d = new TF1("fitfunc", exclusive_xsec_fitfunc, fXmin, fXmax, 1+kNNGFitParams);
     }
  }

  // set limits for fit parameters & fix the others

  LOG("NuVld", pDEBUG) << "Initializing/Fixing NeuGEN fit parameters";

  if(fit_norm) this->InitXSecNormFitParameters();
  else         this->InitFitParameters();

  // get selected data as a TGraph with errors

  TGraphAsymmErrors * gr = 0;

  if( fScaleWithE ) gr = user_data->NuXSec()->GetGraph("all-noE-scale-with-E");
  else gr = user_data->NuXSec()->GetGraph("all-noE");

  // fit graph

  gr->Fit(fFunc1d,"R");

  LOG("NuVld", pDEBUG) << "fitting: ...done";
}
//______________________________________________________________________________
void GuiFitKernel::DoFloatingNormFit(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  // get selected data as a TGraph with errors

  MultiGraph * mgr = user_data->NuXSec()->GetMultiGraph("all-noE-scale-with-E");

  // create a MINUIT fitter with npf params
  // npf = number of fit params (NeuGEN prms + 1 relative norm factor / data-set)
  //       + xsec norm factor (fixed)

  int npf = kNNGFitParams + mgr->NGraphs() + 1;

  TMinuit * minuit = new TMinuit(npf);

  // find out if we are fitting xsecs for an exclusive or inclusive channel

  NeuGenCards * cards = NeuGenCards::Instance();
  bool is_inclusive  = cards->CurrInputs()->Inclusive();

  if( fScaleWithE ) {

     if(is_inclusive) {
         LOG("NuVld", pWARN) << "need to set the fitting for inclusive  channels too";
     }
     else minuit->SetFCN(fcn_mgr_floating_norm_e);

  } else {

     if(is_inclusive) {
         LOG("NuVld", pWARN) << "need to set the fitting for inclusive  channels too";
     }
     else minuit->SetFCN(fcn_mgr_floating_norm);
  }

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist ,1,ierflg);

  this->InitMinuitFitParameters(minuit); // also - need to set the step size...

  // MINUIT's minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  // Printing MINUIT results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit->mnprin(3,amin);

  LOG("NuVld", pDEBUG) << "fitting: ...done";
}
//______________________________________________________________________________
bool GuiFitKernel::MCParamScanning(void)
{
// MC scanning

  if( fNGFP->NFittedParams() <= 0 ) {
     LOG("NuVld", pWARN) << "No parameters to scan";
     return false;
  }

  if(lowb)  delete lowb;  lowb  = 0;
  if(highb) delete highb; highb = 0;

  RandomGen * rnd = RandomGen::Instance();

  NeuGenCards * cards = NeuGenCards::Instance();

  NeuGenConfig * init_config = new NeuGenConfig(cards->CurrConfig()); // orig
  NeuGenWrapper neugen( init_config );

  int           nbins = cards->CurrInputs()->NBins();
  float         emin  = cards->CurrInputs()->EnergyMin();
  float         emax  = cards->CurrInputs()->EnergyMax();
  NGInteraction intr  = cards->CurrInputs()->GetInteraction();
  NGFinalState  fs    = cards->CurrInputs()->GetFinalState();
  NeuGenCuts    cuts  = cards->CurrInputs()->GetCuts();

  bool is_inclusive = cards->CurrInputs()->Inclusive();

  int nmc = 50 * ( fNGFP->NFittedParams() );

  for(int imc = 0; imc < nmc; imc++) {

     LOG("NuVld", pDEBUG) << "MC-iteration: " << imc;

     //-- set a point in the physics parameter space

     for(int ip = 0; ip < kNNGFitParams; ip++) {

        if( fNGFP->IsFitted(ip) ) {

           float min  = (float) fNGFP->RangeMin(ip);
           float max  = (float) fNGFP->RangeMax(ip);

           assert(max>min);

           float x     = rnd->RndGen().Rndm();
           float param = min + x*(max-min);

           LOG("NuVld", pDEBUG) << "Rndm = " << x << ", Param(" << ip << ")=" << param;

           update_neugen_config(param, ip);

        }//is-scanned
     }//ip

     // set current params

     LOG("NuVld", pDEBUG) << "Reconfiguring NeuGEN";

     neugen.Reconfigure( cards->CurrConfig() );

     // get NeuGEN prediction

     Spline * xs_vs_e = 0;

     if(is_inclusive)
        xs_vs_e = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
     else
        xs_vs_e = neugen.ExclusiveXSecSpline(
                                        emin, emax, nbins, &intr, &fs, &cuts);

     // update low/high xsec boundaries

     LOG("NuVld", pDEBUG) << "Updating boundaries";

     if(imc==0) {
       lowb  = new TGraph( * xs_vs_e->GetAsTGraph(1000,fScaleWithE) );
       highb = new TGraph( * xs_vs_e->GetAsTGraph(1000,fScaleWithE) );
     } else {

       for(int ie=0; ie<lowb->GetN(); ie++) {

           double E  = lowb->GetX()[ie];
           double xs = xs_vs_e->Evaluate(E);

           if(fScaleWithE && E>0) xs/=E;

           lowb  -> GetY()[ie] = TMath::Min( xs, lowb  -> GetY()[ie] );
           highb -> GetY()[ie] = TMath::Max( xs, highb -> GetY()[ie] );
       }
    }
  }// imc

  cards->CurrConfig()->Copy(init_config); // restore initial configuration
  delete init_config;

  return true;
}
//______________________________________________________________________________
bool GuiFitKernel::ChisqScan1D(void)
{
  if( fNGFP->NFittedParams() != 1 ) {
     LOG("NuVld", pWARN) << "Number of scanned params != 1, in 1D-Scan!";

     if(chisq1d) delete chisq1d;
     chisq1d = 0;

     return false;
  }

  NeuGenCards * cards = NeuGenCards::Instance();

  NeuGenConfig * init_config = new NeuGenConfig(cards->CurrConfig()); // orig
  NeuGenWrapper neugen( init_config );

  int           nbins = cards->CurrInputs()->NBins();
  float         emin  = cards->CurrInputs()->EnergyMin();
  float         emax  = cards->CurrInputs()->EnergyMax();
  NGInteraction intr  = cards->CurrInputs()->GetInteraction();
  NGFinalState  fs    = cards->CurrInputs()->GetFinalState();
  NeuGenCuts    cuts  = cards->CurrInputs()->GetCuts();

  bool is_inclusive = cards->CurrInputs()->Inclusive();

  Spline * ref = 0, * cur = 0;

  if(is_inclusive)
        ref = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
  else
        ref = neugen.ExclusiveXSecSpline(
                            emin, emax, nbins, &intr, &fs, &cuts);

  double * param_values = 0;
  double * chisq_values = 0;

  bool prm_found = false;

  for(int ip = 0; ip < kNNGFitParams; ip++) {

    if(fNGFP->IsFitted(ip) && !prm_found) {

      float min  = (float) fNGFP->RangeMin(ip);
      float max  = (float) fNGFP->RangeMax(ip);
      float step = (float) fNGFP->Step(ip);

      prm_found = true;

      assert(max>min);

      if(step<=0) {
         LOG("NuVld", pDEBUG) << "step = " << step << "! --> setting default";

         step = 0.04 * (max-min);
      }

      int n = int((max-min)/step);

      LOG("NuVld", pINFO) << "n = " << n;

      assert(n>0);

      param_values = new double[n];
      chisq_values = new double[n];

      for(int i=0; i<=n; i++) {

          // set current param

          float param = min + i * step;

          update_neugen_config(param, ip);

          neugen.Reconfigure( cards->CurrConfig() );

          // get NeuGEN prediction

          if(is_inclusive)
                cur = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
          else
                cur = neugen.ExclusiveXSecSpline(
                                    emin, emax, nbins, &intr, &fs, &cuts);

          // compute chisq

          double chisq = 0;

          float dE = (emax-emin)/(nbins-1);

          for(int ie = 0; ie < nbins; ie++) {

             float E = emin + ie * dE;

             chisq += TMath::Power( ref->Evaluate(E) - cur->Evaluate(E), 2 );
          }//E

          LOG("NuVld", pINFO) << "param = " << param << " - chisq = " << chisq;

          param_values[i] = param;
          chisq_values[i] = chisq;

      } //param values

      if(chisq1d) delete chisq1d;

      chisq1d = new TGraph(n,param_values,chisq_values);

    } //is 1st scanned param
  }//params

  cards->CurrConfig()->Copy(init_config); // restore initial configuration
  delete init_config;

  return true;
}
//______________________________________________________________________________
bool GuiFitKernel::ChisqScan2D(void)
{
  if( fNGFP->NFittedParams() != 2 ) {
     LOG("NuVld", pWARN) << "Number of scanned params != 2, in 2D-Scan!";

     if(chisq2d) delete chisq2d;
     chisq2d = 0;

     return false;
  }

  NeuGenCards * cards = NeuGenCards::Instance();

  NeuGenConfig * init_config = new NeuGenConfig(cards->CurrConfig()); // orig
  NeuGenWrapper neugen( init_config );

  int           nbins = cards->CurrInputs()->NBins();
  float         emin  = cards->CurrInputs()->EnergyMin();
  float         emax  = cards->CurrInputs()->EnergyMax();
  NGInteraction intr  = cards->CurrInputs()->GetInteraction();
  NGFinalState  fs    = cards->CurrInputs()->GetFinalState();
  NeuGenCuts    cuts  = cards->CurrInputs()->GetCuts();

  bool is_inclusive = cards->CurrInputs()->Inclusive();

  Spline * ref = 0, * cur = 0;

  if(is_inclusive)
        ref = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
  else
        ref = neugen.ExclusiveXSecSpline(
                            emin, emax, nbins, &intr, &fs, &cuts);

  int   ip1 = -1, ip2 = -1;
  float min1 = 0, max1 = 0, step1 = 0, min2 = 0, max2 = 0, step2 = 0;

  bool prm_found = false;
  for(int ip = 0; ip < kNNGFitParams; ip++) {

    if(fNGFP->IsFitted(ip) && !prm_found) {

      min1  = (float) fNGFP->RangeMin(ip);
      max1  = (float) fNGFP->RangeMax(ip);
      step1 = (float) fNGFP->Step(ip);

      ip1 = ip;
      prm_found = true;
    }
  }

  prm_found = false;
  for(int ip = 0; ip < kNNGFitParams; ip++) {

    if(fNGFP->IsFitted(ip) && !prm_found && ip!=ip1) {

      min2  = (float) fNGFP->RangeMin(ip);
      max2  = (float) fNGFP->RangeMax(ip);
      step2 = (float) fNGFP->Step(ip);

      ip2 = ip;
      prm_found = true;
    }
  }
  assert(step1>0 && step2>0);
  assert(max1>min1 && max2>min2);

  int n1 = int( (max1-min1)/step1 );
  int n2 = int( (max2-min2)/step2 );

  assert(n1>0);
  assert(n2>0);

  if(chisq2d) delete chisq2d;
  chisq2d = new TH2F("chisq2d","",n1,min1,max1,n2,min2,max2);

  for(int i=0; i<=n1; i++) {
    for(int j=0; j<=n2; j++) {

          // set current param

          float param1 = min1 + i * step1;
          float param2 = min2 + j * step2;

          update_neugen_config(param1, ip1);
          update_neugen_config(param2, ip2);

          neugen.Reconfigure( cards->CurrConfig() );

          // get NeuGEN prediction

          if(is_inclusive)
                cur = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
          else
                cur = neugen.ExclusiveXSecSpline(
                                    emin, emax, nbins, &intr, &fs, &cuts);

          // compute chisq

          double chisq = 0;

          float dE = (emax-emin)/(nbins-1);

          for(int ie = 0; ie < nbins; ie++) {

             float E = emin + ie * dE;

             chisq += TMath::Power( ref->Evaluate(E) - cur->Evaluate(E), 2 );
          }//E

          LOG("NuVld", pINFO) << "params = ("
                          << param1 << ", " << param2 << " - chisq = " << chisq;

          chisq2d->Fill(param1+0.1*step1, param2+0.1*step2, chisq);
    } //param values
  } //param values

  cards->CurrConfig()->Copy(init_config); // restore initial configuration
  delete init_config;

  return true;
}
//______________________________________________________________________________
void GuiFitKernel::InitFitParameters(void)
{
  for(int i = 0; i < kNNGFitParams; i++) {

    fFunc1d->SetParName( i, fNGFP->ParamAsString(i).c_str() );

    if(fNGFP->IsFitted(i))
    {
      // As set in current neugen config
      float value = this->CurrFitParamValue(i);

      // Range as set in neugen fit params dialog
      float min = (float) fNGFP->RangeMin(i);
      float max = (float) fNGFP->RangeMax(i);

      LOG("NuVld", pINFO)
           << "** Setting fit parameter " << i << "(" << fNGFP->ParamAsString(i)
           << ") value = " << value << ", range = [" << min << ", " << max <<"]";

      fFunc1d->SetParameter(i, value);
      fFunc1d->SetParLimits(i, min, max);

    } else {

      // As set in current neugen config
      float value = this->CurrFitParamValue(i);

      LOG("NuVld", pINFO)
           << "** Fixing fit parameter " << i << "(" << fNGFP->ParamAsString(i)
           << ") value := " << value;

      fFunc1d->FixParameter(i, value);
    }
  }//i

  fFunc1d->SetParName(kNNGFitParams, "XSEC NORM");
  fFunc1d->FixParameter(kNNGFitParams, 1.);
}
//______________________________________________________________________________
void GuiFitKernel::InitXSecNormFitParameters(void)
{
  // fix all physics parameters to best values

  for(int i = 0; i < kNNGFitParams; i++) {

    fFunc1d->SetParName( i, fNGFP->ParamAsString(i).c_str() );

    // As set in current neugen config
    float value = this->CurrFitParamValue(i);

    LOG("NuVld", pINFO)
           << "** Fixing fit parameter " << i << "(" << fNGFP->ParamAsString(i)
           << ") value := " << value;

    fFunc1d->FixParameter(i, value);
  }

  // allow an overall normalization parameter

  fFunc1d->SetParName(kNNGFitParams, "XSEC NORM");
  fFunc1d->SetParameter(kNNGFitParams, 1.);
  fFunc1d->SetParLimits(kNNGFitParams, 0.,2.);
}
//______________________________________________________________________________
void GuiFitKernel::InitMinuitFitParameters(TMinuit * minuit)
{
  int ierrflag = 0;

  NuVldUserData * user_data = NuVldUserData::Instance();

  MultiGraph * mgr = user_data->NuXSec()->GetMultiGraph("all-noE-scale-with-E");

  int npf = kNNGFitParams + mgr->NGraphs() + 1;

  for(int i = 0; i < kNNGFitParams; i++) {

    float value = this->CurrFitParamValue(i);

    // Range as set in neugen fit params dialog
    float min  = (float) fNGFP->RangeMin(i);
    float max  = (float) fNGFP->RangeMax(i);
    float step = (float) fNGFP->Step(i);

    if(step < 0) step = 0.01;

    LOG("NuVld", pINFO)
           << "** Setting fit parameter " << i << "(" << fNGFP->ParamAsString(i)
           << ") value = " << value << ", range = [" << min << ", " << max <<"]";

    minuit->mnparm(i, fNGFP->ParamAsString(i).c_str(), value, step, min,max, ierrflag);
  }

  fFunc1d->SetParName(kNNGFitParams, "XSEC NORM");
  fFunc1d->FixParameter(kNNGFitParams, 1.);

  int igraph = 0;
  for(int i = kNNGFitParams+1; i < npf; i++) {
    minuit->mnparm(i, mgr->GetLegendEntry(igraph++), 1.0, 0.03, 0.8, 1.2,ierrflag);
  }

  for(int i = 0; i < kNNGFitParams; i++) {
    if( ! fNGFP->IsFitted(i)) {
      LOG("NuVld", pINFO)
           << "** Fixing fit parameter " << i << "(" << fNGFP->ParamAsString(i);
      minuit->FixParameter(i);
    }
  }
}
//______________________________________________________________________________
float GuiFitKernel::CurrFitParamValue(int iparam)
{
  NeuGenCards * cards = NeuGenCards::Instance();

  switch(iparam) {

  case ( kNgfMaQel       ): return cards->CurrConfig()->MaQel();          break;
  case ( kNgfMaRes       ): return cards->CurrConfig()->MaRes();          break;
  case ( kNgfMaCoh       ): return cards->CurrConfig()->MaCoh();          break;
  case ( kNgfQelFa0      ): return cards->CurrConfig()->Fa0Qel();         break;
  case ( kNgfQelEta      ): return cards->CurrConfig()->EtaQel();         break;
  case ( kNgfResOmega    ): return cards->CurrConfig()->OmegaRes();       break;
  case ( kNgfResZ        ): return cards->CurrConfig()->ZRes();           break;
  case ( kNgfCohR0       ): return cards->CurrConfig()->R0Coh();          break;
  case ( kNgfCohREI      ): return cards->CurrConfig()->REICoh();         break;
  case ( kNgfKnoBvp      ): return cards->CurrConfig()->KnoB(e_vp);       break;
  case ( kNgfKnoBvn      ): return cards->CurrConfig()->KnoB(e_vn);       break;
  case ( kNgfKnoBvbp     ): return cards->CurrConfig()->KnoB(e_vbp);      break;
  case ( kNgfKnoBvbn     ): return cards->CurrConfig()->KnoB(e_vbn);      break;
  case ( kNgfKnoAvp      ): return cards->CurrConfig()->KnoA(e_vp);       break;
  case ( kNgfKnoAvn      ): return cards->CurrConfig()->KnoA(e_vn);       break;
  case ( kNgfKnoAvbp     ): return cards->CurrConfig()->KnoA(e_vbp);      break;
  case ( kNgfKnoAvbn     ): return cards->CurrConfig()->KnoA(e_vbn);      break;
  case ( kNgfDisResM2vp  ): return cards->CurrConfig()->DisRes(2, e_vp);  break;
  case ( kNgfDisResM3vp  ): return cards->CurrConfig()->DisRes(3, e_vp);  break;
  case ( kNgfDisResM2vn  ): return cards->CurrConfig()->DisRes(2, e_vn);  break;
  case ( kNgfDisResM3vn  ): return cards->CurrConfig()->DisRes(3, e_vn);  break;
  case ( kNgfDisResM2vbp ): return cards->CurrConfig()->DisRes(2, e_vbp); break;
  case ( kNgfDisResM3vbp ): return cards->CurrConfig()->DisRes(3, e_vbp); break;
  case ( kNgfDisResM2vbn ): return cards->CurrConfig()->DisRes(2, e_vbn); break;
  case ( kNgfDisResM3vbn ): return cards->CurrConfig()->DisRes(3, e_vbn); break;

  default :                 return 0;
  }
  return 0;
}
//______________________________________________________________________________
TGraph * GuiFitKernel::GetResidualsAsGraph(void)
{
  if(!fFunc1d) return 0;

  //-- getting the fitted data as TGraphAsymmErrors

  NuVldUserData * user_data = NuVldUserData::Instance();

  TGraphAsymmErrors * gr = user_data->NuXSec()->GetGraph("all-noE");

  if(!gr) return 0;

  //-- computing the residuals

  double * x = new double[ gr->GetN() ];
  double * y = new double[ gr->GetN() ];

  for(int i = 0; i < gr->GetN() ; i++) {

     x[i] = gr->GetX()[i];
     y[i] = gr->GetY()[i] - fFunc1d->Eval( x[i] );
  }

  TGraph * residuals = new TGraph( gr->GetN(), x, y );

  delete [] x;
  delete [] y;

  return residuals;
}
//______________________________________________________________________________
*/
