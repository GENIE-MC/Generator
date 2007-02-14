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

#include <cstdlib>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/FGMBodekRitchie.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
FGMBodekRitchie::FGMBodekRitchie() :
NuclearModelI("genie::FGMBodekRitchie")
{

}
//____________________________________________________________________________
FGMBodekRitchie::FGMBodekRitchie(string config) :
NuclearModelI("genie::FGMBodekRitchie", config)
{

}
//____________________________________________________________________________
FGMBodekRitchie::~FGMBodekRitchie()
{

}
//____________________________________________________________________________
bool FGMBodekRitchie::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  TH1D * prob = this->ProbDistro(target);
  if(!prob) {
    LOG("BodekRitchie", pFATAL)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }

  double p = prob->GetRandom();
  LOG("BodekRitchie", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;  

  fCurrRemovalEnergy = fNucRmvE;;
  fCurrMomentum.SetXYZ(px,py,pz);

  return true;
}
//____________________________________________________________________________
double FGMBodekRitchie::Prob(double p, double w, const Target & target) const
{
  if(w<0) {
     TH1D * prob = this->ProbDistro(target);
     int bin = prob->FindBin(p);
     double y  = prob->GetBinContent(bin);
     double dx = prob->GetBinWidth(bin);
     double p  = y*dx;
     return p;
  }
  return 1;
}
//____________________________________________________________________________
TH1D * FGMBodekRitchie::ProbDistro(const Target & target) const
{
  //-- return stored /if already computed/
  map<string, TH1D*>::iterator it = fProbDistroMap.find(target.AsString());
  if(it != fProbDistroMap.end()) return it->second;

  LOG("BodekRitchie", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("BodekRitchie", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.HitNucPdg();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  double Z = (double) target.Z();
  double N = (double) target.N();
  double A = (double) target.A();

  //-- look up the Fermi momentum for this Target
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable(fKFTable);
  double KF = kft->FindClosestKF(target_pdgc, nucleon_pdgc);
  LOG("BodekRitchie", pNOTICE) << "KF = " << KF;

  //-- correct the Fermi momentum for the struck nucleon
  assert(target.HitNucIsSet());
  bool is_p = pdg::IsProton(nucleon_pdgc);
  if(is_p) KF *= TMath::Power( 2*Z/A, 1./3.);
  else     KF *= TMath::Power( 2*N/A, 1./3.);
  LOG("BodekRitchie", pINFO) << "Corrected KF = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;
  double R  = 1. / (1.- KF/4.);
  LOG("BodekRitchie", pDEBUG) << "a  = " << a;
  LOG("BodekRitchie", pDEBUG) << "C  = " << C;
  LOG("BodekRitchie", pDEBUG) << "R  = " << R;

  //-- create the probability distribution

  TH1D * prob = new TH1D("", "", fNPBins, 0, fPMax);

  double dp = fPMax / (fNPBins-1);
  double iC = (C>0) ? 1./C : 0.;
  double kfa_pi_2 = TMath::Power(KF*a/kPi,2);

  for(int i = 0; i < fNPBins; i++) {
     double p  = i * dp;
     double p2 = TMath::Power(p,2);

     // calculate |phi(p)|^2
     double phi2 = 0;
     if (p <= KF)
        phi2 = iC * (1. - 6.*kfa_pi_2);
     else if ( p > KF && p < fPCutOff)
        phi2 = iC * (2*R*kfa_pi_2*TMath::Power(KF/p,4.));

     // calculate probability density : dProbability/dp
     double dP_dp = 4*kPi * p2 * phi2;
     LOG("BodekRitchie", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  //-- store
  fProbDistroMap.insert(
      map<string, TH1D*>::value_type(target.AsString(),prob));

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
void FGMBodekRitchie::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FGMBodekRitchie::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void FGMBodekRitchie::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fKFTable = fConfig->GetStringDef ("FermiMomentumTable", 
                                    gc->GetString("FermiMomentumTable"));
  fNPBins  = fConfig->GetIntDef    ("Momentum-NumBins",  400);
  fPMax    = fConfig->GetDoubleDef ("Momentum-Max",     10.0);
  fPCutOff = fConfig->GetDoubleDef ("Momentum-CutOff",  
                                    gc->GetDouble("RFG-Momentum-CutOff"));
  fNucRmvE = fConfig->GetDoubleDef ("NucRemovalE",  
                                    gc->GetDouble("RFG-NucRemovalE"));

  fNucRmvE = TMath::Max(fNucRmvE, 0.);

  assert(fNPBins > 1 && fPMax > 0 && fPCutOff > 0 && fPCutOff < fPMax);
}
//____________________________________________________________________________

