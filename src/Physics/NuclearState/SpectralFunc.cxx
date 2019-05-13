//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 01, 2012 - CA
   Pick spectral function data from $GENIE/data/evgen/nucl/spectral_functions
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
SpectralFunc::SpectralFunc() :
NuclearModelI("genie::SpectralFunc")
{
  fSfFe56 = 0;
  fSfC12  = 0;
}
//____________________________________________________________________________
SpectralFunc::SpectralFunc(string config) :
NuclearModelI("genie::SpectralFunc", config)
{
  fSfFe56 = 0;
  fSfC12  = 0;
}
//____________________________________________________________________________
SpectralFunc::~SpectralFunc()
{
//if (fSfFe56) delete fSfFe56;
//if (fSfC12 ) delete fSfC12;
}
//____________________________________________________________________________
bool SpectralFunc::GenerateNucleon(const Target & target) const
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

  LOG("SpectralFunc", pINFO) << "Momentum range = ["   << kmin << ", " << kmax << "]"; 
  LOG("SpectralFunc", pINFO) << "Rmv energy range = [" << wmin << ", " << wmax << "]";

  RandomGen * rnd = RandomGen::Instance();

  unsigned int niter = 0;
  while(1) {
    if(niter > kRjMaxIterations) {
       LOG("SpectralFunc", pWARN) 
           << "Couldn't generate a hit nucleon after " << niter << " iterations";
       return false;
    }
    niter++;

    // random pair
    double kc = kmin + dk * rnd->RndGen().Rndm();
    double wc = wmin + dw * rnd->RndGen().Rndm();
    LOG("SpectralFunc", pINFO) << "Trying p = " << kc << ", w = " << wc;

    // accept/reject
    double prob  = this->Prob(kc,wc, target);
    double probg = probmax * rnd->RndGen().Rndm();
    bool accept = (probg < prob);
    if(!accept) continue;

    LOG("SpectralFunc", pINFO) << "|p,nucleon| = " << kc; 
    LOG("SpectralFunc", pINFO) << "|w,nucleon| = " << wc;

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
double SpectralFunc::Prob(
                         double p, double w, const Target & target) const
{
  TGraph2D * sf = this->SelectSpectralFunction(target);
  if(!sf) return 0;

  return sf->Interpolate(p,w);
}
//____________________________________________________________________________
void SpectralFunc::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::LoadConfig(void)
{
  LOG("SpectralFunc", pDEBUG) << "Loading Benhar et al. spectral functions";

  string data_dir =
        string(gSystem->Getenv("GENIE")) + 
        string("/data/evgen/nucl/spectral_functions/");
  string c12file  = data_dir + "benhar-sf-12c.data";
  string fe56file = data_dir + "benhar-sf-56fe.data";

  TNtupleD sfdata_fe56("sfdata_fe56","","k:e:prob");
  TNtupleD sfdata_c12 ("sfdata_c12", "","k:e:prob");

  sfdata_fe56.ReadFile ( fe56file.c_str() );
  sfdata_c12. ReadFile ( c12file.c_str () );

  LOG("SpectralFunc", pDEBUG) << "Loaded " << sfdata_fe56.GetEntries() << " Fe56 points";
  LOG("SpectralFunc", pDEBUG) << "Loaded " << sfdata_c12.GetEntries()  << " C12 points";

  if (fSfFe56) delete fSfFe56;
  if (fSfC12 ) delete fSfC12;

  fSfFe56 = this->Convert2Graph(sfdata_fe56);
  fSfC12  = this->Convert2Graph(sfdata_c12);

  fSfFe56->SetName("sf_fe56");
  fSfC12 ->SetName("sf_c12");
}
//____________________________________________________________________________
TGraph2D * SpectralFunc::Convert2Graph(TNtupleD & sfdata) const
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
TGraph2D * SpectralFunc::SelectSpectralFunction(const Target & t) const
{
  TGraph2D * sf = 0;
  int pdgc = t.Pdg();

  if      (pdgc == kPdgTgtC12)  sf = fSfC12;
  else if (pdgc == kPdgTgtFe56) sf = fSfFe56;
  else {
    LOG("SpectralFunc", pERROR) 
     << "** The spectral function for target " << pdgc << " isn't available";
  }
  if(!sf) {
    LOG("SpectralFunc", pERROR) << "** Null spectral function";
  }
  return sf;
}
//____________________________________________________________________________
