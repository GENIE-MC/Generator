//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - June 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TSystem.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Messenger/Messenger.h"
#include "Nuclear/BenharSpectralFunc1D.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
BenharSpectralFunc1D::BenharSpectralFunc1D() :
NuclearModelI("genie::BenharSpectralFunc1D")
{
  fSfC12_k     = 0;   
  fSfFe56_k    = 0; 
  fMaxC12Prob  = 0;
  fMaxFe56Prob = 0;
}
//____________________________________________________________________________
BenharSpectralFunc1D::BenharSpectralFunc1D(string config) :
NuclearModelI("genie::BenharSpectralFunc1D", config)
{
  fSfC12_k     = 0;   
  fSfFe56_k    = 0; 
  fMaxC12Prob  = 0;
  fMaxFe56Prob = 0;
}
//____________________________________________________________________________
BenharSpectralFunc1D::~BenharSpectralFunc1D()
{
  if (fSfC12_k ) delete fSfC12_k;
  if (fSfFe56_k) delete fSfFe56_k;
}
//____________________________________________________________________________
bool BenharSpectralFunc1D::GenerateNucleon(const Target & target) const
{
  Spline * distrib = this->SelectMomentumDistrib(target);
  if(!distrib) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0.,0.,0.);
    return false;
  }

  RandomGen * rnd = RandomGen::Instance();

  double max = this->MaxProb(target);
  LOG("BenharSF", pDEBUG) << "Max probability = " << max;

  double p = 0;
  unsigned int niter = 0;
  while(1) {
    if(niter > kRjMaxIterations) {
       LOG("BenharSF", pWARN)
        << "Couldn't generate a hit nucleon after " << niter << " iterations";
       return false;
    }
    niter++;

    p = rnd->RndGen().Rndm();
    double prob  = distrib->Evaluate(p);
    double probg = max * rnd->RndGen().Rndm();
    LOG("BenharSF", pDEBUG) << "Trying p = " << p << " / prob = " << prob;

    bool accept = (probg<prob);
    if(accept) break;
  }

  LOG("BenharSF", pINFO) << "|p,nucleon| = " << p;

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  fCurrRemovalEnergy = 0.025;
  fCurrMomentum.SetXYZ(px,py,pz);

  return true;
}
//____________________________________________________________________________
double BenharSpectralFunc1D::Prob(
                  double p, double /*w*/, const Target & target) const
{
  Spline * distrib = this->SelectMomentumDistrib(target);
  if(!distrib) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0.,0.,0.);
    return false;
  }
  double prob  = distrib->Evaluate(p*1000);
  return prob;
}
//____________________________________________________________________________
void BenharSpectralFunc1D::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BenharSpectralFunc1D::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void BenharSpectralFunc1D::LoadConfig(void)
{
  LOG("BenharSF", pDEBUG) << "Loading coonfiguration for BenharSpectralFunc1D";

  if (fSfC12_k ) delete fSfC12_k;
  if (fSfFe56_k) delete fSfFe56_k;
  
  // load 
  string data_dir = 
        string(gSystem->Getenv("GENIE")) + string("/data/spectral_functions/");
  string c12file  = data_dir + "benhar-sf1dk-12c.data";
  string fe56file = data_dir + "benhar-sf1dk-56fe.data";

  fSfC12_k  = new Spline( c12file  );
  fSfFe56_k = new Spline( fe56file );

  // scan for max.
  fMaxC12Prob  = 0;
  fMaxFe56Prob = 0;
  int    n    = 100;
  double pmin =  0.000;
  double dp   =  0.010;
  for(int i=0; i<n; i++) {
    double p = pmin + i*dp;
    fMaxC12Prob  = TMath::Max(fMaxC12Prob,  fSfC12_k ->Evaluate(p));
    fMaxFe56Prob = TMath::Max(fMaxFe56Prob, fSfFe56_k->Evaluate(p));
  }
}
//____________________________________________________________________________
Spline * BenharSpectralFunc1D::SelectMomentumDistrib(const Target & t) const
{
  Spline * distrib = 0;
  int pdgc = t.Pdg();
    
  if      (pdgc == kPdgTgtC12)  distrib = fSfC12_k;
  else if (pdgc == kPdgTgtFe56) distrib = fSfFe56_k;
  else {
    LOG("BenharSF", pERROR)   
    << "** The momentum distribution for target " << pdgc << " isn't available";
  }
  if(!distrib) {
    LOG("BenharSF", pERROR) << "** Null momentum distribution";
  }                      
  return distrib;
}
//____________________________________________________________________________
double BenharSpectralFunc1D::MaxProb(const Target & t) const
{
  int pdgc = t.Pdg();
    
  if      (pdgc == kPdgTgtC12)  return fMaxC12Prob;
  else if (pdgc == kPdgTgtFe56) return fMaxFe56Prob;
  else                          return 0;
}
//____________________________________________________________________________


