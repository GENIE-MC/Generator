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

#include <TFile.h>

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/BenharSpectralFunc1D.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BenharSpectralFunc1D::BenharSpectralFunc1D() :
NuclearModelI("genie::BenharSpectralFunc1D")
{

}
//____________________________________________________________________________
BenharSpectralFunc1D::BenharSpectralFunc1D(string config) :
NuclearModelI("genie::BenharSpectralFunc1D", config)
{

}
//____________________________________________________________________________
BenharSpectralFunc1D::~BenharSpectralFunc1D()
{

}
//____________________________________________________________________________
bool BenharSpectralFunc1D::GenerateNucleon(const Target & tgt) const
{
  double p = fHisto->GetRandom();
  LOG("BenharSF1D", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

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
double BenharSpectralFunc1D::Prob(double p, double w, const Target & target) const
{
  return 1;
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
  LOG("BenharSF1D", pNOTICE)
             << "Loading coonfiguration for DytmanSF";

  // NOTE: temporary
  //
  TFile f("/home/costas/mywork/GENIE/src/Nuclear/sf56.root","read");

  TH1D * h = (TH1D*)f.Get("sfFe56_k");
  LOG("BenharSF1D", pNOTICE) << "mean = " << h->GetMean();

  fHisto = new TH1D(*h);
  fHisto->SetDirectory(0);
  fHisto->Scale( 1.0 / fHisto->Integral("width") );

  f.Close(); 

  LOG("BenharSF1D", pNOTICE) << "mean = " << fHisto->GetMean();
}
//____________________________________________________________________________

