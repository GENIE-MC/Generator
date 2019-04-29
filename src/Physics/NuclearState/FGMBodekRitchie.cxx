//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2008 - CA
   Call SetDirectory(0) at the temp momentum distribution histogram to stop
   it from being automatically written out at the event file.
 @ Jun 18, 2008 - CA
   Deallocate the momentum distribution histograms map at dtor
*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FGMBodekRitchie.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

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
  map<string, TH1D*>::iterator iter = fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.end(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//____________________________________________________________________________
bool FGMBodekRitchie::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  //-- set fermi momentum vector
  //
  TH1D * prob = this->ProbDistro(target);
  if ( ! prob ) {
    LOG("BodekRitchie", pNOTICE)
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

  fCurrMomentum.SetXYZ(px,py,pz);

  //-- set removal energy
  //
  if ( target.A() < 6 || ! fUseParametrization )
  {
     int Z = target.Z();
     map<int,double>::const_iterator it = fNucRmvE.find(Z);
     if(it != fNucRmvE.end()) fCurrRemovalEnergy = it->second;
     else fCurrRemovalEnergy = nuclear::BindEnergyPerNucleon(target);
  }
  else {
     fCurrRemovalEnergy = nuclear::BindEnergyPerNucleonParametrization(target);
  }

  return true;
}
//____________________________________________________________________________
double FGMBodekRitchie::Prob(double mom, double w, const Target & target) const
{
  if(w<0) {
     TH1D * prob_distr = this->ProbDistro(target);
     int bin = prob_distr->FindBin(mom);
     double y  = prob_distr->GetBinContent(bin);
     double dx = prob_distr->GetBinWidth(bin);
     double prob = y*dx;
     return prob;
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
  double KF;
  if (A<6 || !fUseParametrization)
  {
     //-- look up the Fermi momentum for this Target
     FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
     const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
     KF = kft->FindClosestKF(target_pdgc, nucleon_pdgc);
     LOG("BodekRitchie", pNOTICE) << "KF = " << KF;
  }
  else
  {
      //-- define the Fermi momentum for this Target
      //
      KF = nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
      //-- correct the Fermi momentum for the struck nucleon
      assert(target.HitNucIsSet());
      bool is_p = pdg::IsProton(nucleon_pdgc);
      if(is_p) KF *= TMath::Power( 2*Z/A, 1./3.);
      else
         KF *= TMath::Power( 2*N/A, 1./3.);
      LOG("BodekRitchie", pINFO) << "Corrected KF = " << KF;
   }




  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;
  double R  = 1. / (1.- KF/fPMax);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekRitchie", pDEBUG) << "a  = " << a;
  LOG("BodekRitchie", pDEBUG) << "C  = " << C;
  LOG("BodekRitchie", pDEBUG) << "R  = " << R;
#endif

  //-- create the probability distribution

  int npbins = (int) (1000*fPMax);
  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);
  prob->SetDirectory(0);

  double dp = fPMax / (npbins-1);
  double iC = (C>0) ? 1./C : 0.;
  double kfa_pi_2 = TMath::Power(KF*a/kPi,2);

  for(int i = 0; i < npbins; i++) {
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
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("BodekRitchie", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  //-- store
  fProbDistroMap.insert(
      map<string, TH1D*>::value_type(target.AsString(),prob));

  return prob;
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
  this->GetParam( "FermiMomentumTable", fKFTable);

  
  // Default value 4.0 from original paper by A. Bodek and J. L. Ritchie. Phys. Rev. D 23, 1070
  this->GetParamDef("MomentumMax", fPMax, 4.0);
  this->GetParam("RFG-MomentumCutOff", fPCutOff);
  this->GetParam("RFG-UseParametrization", fUseParametrization);

  assert(fPMax > 0 && fPCutOff > 0 && fPCutOff < fPMax);

  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
  const std::string keyStart = "RFG-NucRemovalE@Pdg=";

  RgIMap entries = GetConfig().GetItemMap();

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string& key = it->first;
    int pdg = 0;
    int Z = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str())) {
      pdg = atoi(key.c_str() + keyStart.size());
      Z = pdg::IonPdgCodeToZ(pdg);
    }
    if (0 != pdg && 0 != Z) {
      ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double eb ;
      GetParam( rgkey, eb ) ;
      eb = TMath::Max(eb, 0.);
      LOG("BodekRitchie", pINFO)
        << "Nucleus: " << pdg << " -> using Eb =  " << eb << " GeV";
      fNucRmvE.insert(map<int,double>::value_type(Z,eb));
    }
  }

  LOG("BodekRitchie", pDEBUG)
    << "Finished LoadConfig";
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  for (map<int,double>::iterator it = fNucRmvE.begin(); it != fNucRmvE.end(); ++it) {
    LOG("BodekRitchie", pDEBUG)
      << "Z = " << (*it).first << "; eb = " << (*it).second;
  }
#endif
}
//____________________________________________________________________________
