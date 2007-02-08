//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - June 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Nuclear/BenharSpectralFunc.h"
#include "PDG/PDGCodes.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
BenharSpectralFunc::BenharSpectralFunc() :
NuclearModelI("genie::BenharSpectralFunc")
{
  fSfFe56 = 0;
  fSfC12  = 0;
}
//____________________________________________________________________________
BenharSpectralFunc::BenharSpectralFunc(string config) :
NuclearModelI("genie::BenharSpectralFunc", config)
{
  fSfFe56 = 0;
  fSfC12  = 0;
}
//____________________________________________________________________________
BenharSpectralFunc::~BenharSpectralFunc()
{
  if (fSfFe56) delete fSfFe56;
  if (fSfC12 ) delete fSfC12;
}
//____________________________________________________________________________
bool BenharSpectralFunc::GenerateNucleon(const Target & target) const
{
  TGraph2D * sf = this->SelectSpectralFunction(target);

  if(!sf) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0.,0.,0.);
    return false;
  }

  double kmin    = sf->GetXmin(); // momentum range
  double kmax    = sf->GetXmax();
  double wmin    = sf->GetYmin(); // removal energy range
  double wmax    = sf->GetYmax();
  double probmax = sf->GetZmax(); // maximum probability
  probmax *= 1.1;

  double dk = kmax - kmin;
  double dw = wmax - wmin;

  LOG("BenharSF", pINFO) << "Momentum range = ["   << kmin << ", " << kmax << "]"; 
  LOG("BenharSF", pINFO) << "Rmv energy range = [" << wmin << ", " << wmax << "]";

  RandomGen * rnd = RandomGen::Instance();

  unsigned int niter = 0;
  while(1) {
    if(niter > kRjMaxIterations) {
       LOG("BenharSF", pWARN) 
           << "Couldn't generate a hit nucleon after " << niter << " iterations";
       return false;
    }
    niter++;

    // random pair
    double kc = kmin + dk * rnd->RndGen().Rndm();
    double wc = wmin + dw * rnd->RndGen().Rndm();
    LOG("BenharSF", pINFO) << "Trying p = " << kc << ", w = " << wc;

    // accept/reject
    double prob  = this->Prob(kc,wc, target);
    double probg = probmax * rnd->RndGen().Rndm();
    bool accept = (probg < prob);
    if(!accept) continue;

    LOG("BenharSF", pINFO) << "|p,nucleon| = " << kc; 
    LOG("BenharSF", pINFO) << "|w,nucleon| = " << wc;

    // generate momentum components
    double costheta = -1. + 2. * rnd->RndGen().Rndm();
    double sintheta = TMath::Sqrt(1.-costheta*costheta);
    double fi       = 2 * kPi * rnd->RndGen().Rndm();
    double cosfi    = TMath::Cos(fi);
    double sinfi    = TMath::Sin(fi);

    double kx = kc*sintheta*cosfi;
    double ky = kc*sintheta*sinfi;
    double kz = kc*costheta;

    // set generated values
    fCurrRemovalEnergy = wc;
    fCurrMomentum.SetXYZ(kx,ky,kz);

    return true;
  }
  return false;
}
//____________________________________________________________________________
double BenharSpectralFunc::Prob(
                         double p, double w, const Target & target) const
{
  TGraph2D * sf = this->SelectSpectralFunction(target);
  if(!sf) return 0;

  return sf->Interpolate(p,w);
}
//____________________________________________________________________________
void BenharSpectralFunc::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BenharSpectralFunc::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void BenharSpectralFunc::LoadConfig(void)
{
  LOG("BenharSF", pDEBUG) << "Loading Benhar et al. spectral functions";

  string data_dir =
        string(gSystem->Getenv("GENIE")) + string("/data/spectral_functions/");
  string c12file  = data_dir + "benhar-sf-12c.data";
  string fe56file = data_dir + "benhar-sf-56fe.data";

  TNtupleD sfdata_fe56("sfdata_fe56","","k:e:prob");
  TNtupleD sfdata_c12 ("sfdata_c12", "","k:e:prob");

  sfdata_fe56.ReadFile ( fe56file.c_str() );
  sfdata_c12. ReadFile ( c12file.c_str () );

  LOG("BenharSF", pDEBUG) << "Loaded " << sfdata_fe56.GetEntries() << " Fe56 points";
  LOG("BenharSF", pDEBUG) << "Loaded " << sfdata_c12.GetEntries()  << " C12 points";

  if (fSfFe56) delete fSfFe56;
  if (fSfC12 ) delete fSfC12;

  fSfFe56 = this->Convert2Graph(sfdata_fe56);
  fSfC12  = this->Convert2Graph(sfdata_c12);

  fSfFe56->SetName("sf_fe56");
  fSfC12 ->SetName("sf_c12");
}
//____________________________________________________________________________
TGraph2D * BenharSpectralFunc::Convert2Graph(TNtupleD & sfdata) const
{
  int np = sfdata.GetEntries();
  TGraph2D * sfgraph = new TGraph2D(np);

  sfdata.Draw("k:e:prob","","GOFF");
  assert(np==sfdata.GetSelectedRows());
  double * k = sfdata.GetV1();
  double * e = sfdata.GetV2();
  double * p = sfdata.GetV3();

  for(int i=0; i<np; i++) {
    double ki = k[i] * (units::MeV/units::GeV); // momentum
    double ei = e[i] * (units::MeV/units::GeV); // removal energy
    double pi = p[i] * TMath::Power(ki,2);      // probabillity
    sfgraph->SetPoint(i, ki,ei,pi);
  }
  return sfgraph;
}
//____________________________________________________________________________
TGraph2D * BenharSpectralFunc::SelectSpectralFunction(const Target & t) const
{
  TGraph2D * sf = 0;
  int pdgc = t.Pdg();

  if      (pdgc == kPdgTgtC12)  sf = fSfC12;
  else if (pdgc == kPdgTgtFe56) sf = fSfFe56;
  else {
    LOG("BenharSF", pERROR) 
     << "** The spectral function for target " << pdgc << " isn't available";
  }
  if(!sf) {
    LOG("BenharSF", pERROR) << "** Null spectral function";
  }
  return sf;
}
//____________________________________________________________________________
