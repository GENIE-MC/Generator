//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Brian Coopersmith, University of Rochester

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/EffectiveSF.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/ConfigIsotopeMapUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace genie::utils::config;

//____________________________________________________________________________
EffectiveSF::EffectiveSF() :
NuclearModelI("genie::EffectiveSF")
{

}
//____________________________________________________________________________
EffectiveSF::EffectiveSF(string config) :
NuclearModelI("genie::EffectiveSF", config)
{

}
//____________________________________________________________________________
// Cleans up memory from probability distribution maps
//____________________________________________________________________________
EffectiveSF::~EffectiveSF()
{
  map<string, TH1D*>::iterator iter = fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.begin(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//____________________________________________________________________________
// Set the removal energy, 3 momentum, and FermiMover interaction type
//____________________________________________________________________________
bool EffectiveSF::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());
  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  //-- set fermi momentum vector
  //

  if ( target.A() > 1 ) {
    TH1D * prob = this->ProbDistro(target);
    if(!prob) {
      LOG("EffectiveSF", pNOTICE)
              << "Null nucleon momentum probability distribution";
      exit(1);
    }

    double p = prob->GetRandom();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "|p,nucleon| = " << p;
#endif

    RandomGen * rnd = RandomGen::Instance();

    double costheta = -1. + 2. * rnd->RndGen().Rndm();
    double sintheta = TMath::Sqrt(1.-costheta*costheta);
    double fi       = 2 * kPi * rnd->RndGen().Rndm();
    double cosfi    = TMath::Cos(fi);
    double sinfi    = TMath::Sin(fi);

    double px = p*sintheta*cosfi;
    double py = p*sintheta*sinfi;
    double pz = p*costheta;

    fCurrMomentum.SetXYZ(px, py, pz);

  }

  //-- set removal energy
  //

  fCurrRemovalEnergy = this->ReturnBindingEnergy(target);
  double f1p1h = this->Returnf1p1h(target);
  // Since TE increases the QE peak via a 2p2h process, we decrease f1p1h
  // in order to increase the 2p2h interaction to account for this enhancement.
  f1p1h /= this->GetTransEnh1p1hMod(target);
  if ( RandomGen::Instance() -> RndGen().Rndm() < f1p1h) {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF1p1h;
  } else if (fEjectSecondNucleon2p2h) {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF2p2h_eject;
  } else {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF2p2h_noeject;
  }

  return true;
}
//____________________________________________________________________________
// Returns the probability of the bin with given momentum. I don't know what w
// is supposed to be, but I copied its implementation from Bodek-Ritchie.
// Implements the interface.
//____________________________________________________________________________
double EffectiveSF::Prob(double mom, double w, const Target & target) const
{
  if(w < 0) {
     TH1D * prob_distr = this->ProbDistro(target);
     int bin = prob_distr->FindBin(mom);
     double y  = prob_distr->GetBinContent(bin);
     double dx = prob_distr->GetBinWidth(bin);
     double prob  = y * dx;
     return prob;
  }
  return 1;
}
//____________________________________________________________________________
// Check the map of nucleons to see if we have a probability distribution to
// compute with.  If not, make one.
//____________________________________________________________________________
TH1D * EffectiveSF::ProbDistro(const Target & target) const
{
  //-- return stored /if already computed/
  map<string, TH1D*>::iterator it = fProbDistroMap.find(target.AsString());
  if(it != fProbDistroMap.end()) return it->second;

  LOG("EffectiveSF", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("EffectiveSF", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  return this->MakeEffectiveSF(target);

}
//____________________________________________________________________________
// If transverse enhancement form factor modification is enabled, we must
// increase the 2p2h contribution to account for the QE peak enhancement.
// This gets that factor based on the target.
//____________________________________________________________________________
double EffectiveSF::GetTransEnh1p1hMod(const Target& target) const {
  double transEnhMod;
	if(GetValueFromNuclearMaps(target, fTransEnh1p1hMods,
	                           fRangeTransEnh1p1hMods, &transEnhMod) &&
	   transEnhMod > 0) {
	  return transEnhMod;
	}
  // If none specified, assume no enhancement to quasielastic peak
  return 1.0;
}
//____________________________________________________________________________
// Makes a momentum distribtuion for the given target using parameters
// from the config file.
//____________________________________________________________________________
TH1D * EffectiveSF::MakeEffectiveSF(const Target & target) const
{
  // First check for individually specified nuclei
  int pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  map<int,vector<double> >::const_iterator it = fProbDistParams.find(pdgc);
  if(it != fProbDistParams.end()) {
    vector<double> v = it->second;
    return this->MakeEffectiveSF(v[0], v[1], v[2], v[3],
                                 v[4], v[5], v[6], target);
  }

  // Then check in the ranges of A
  map<pair<int, int>, vector<double> >::const_iterator range_it = fRangeProbDistParams.begin();
  for(; range_it != fRangeProbDistParams.end(); ++range_it) {
    if (target.A() >= range_it->first.first && target.A() <= range_it->first.second) {
      vector<double> v = range_it->second;
      return this->MakeEffectiveSF(v[0], v[1], v[2], v[3],
                                   v[4], v[5], v[6], target);
    }
  }

  return NULL;
}
//____________________________________________________________________________
// Makes a momentum distribution using the factors below (see reference) and
// inserts it into the nucleus/momentum distribution map.
//____________________________________________________________________________
TH1D * EffectiveSF::MakeEffectiveSF(double bs, double bp, double alpha,
                                    double beta, double c1, double c2,
                                    double c3, const Target & target) const
{
  //-- create the probability distribution
  int npbins = (int) (1000 * fPMax);

  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);
  prob->SetDirectory(0);

  double dp = fPMax / (npbins-1);
  for(int i = 0; i < npbins; i++) {
    double p  = i * dp;
    double y = p / 0.197;
    double as = c1 * exp(-pow(bs*y,2));
    double ap = c2 * pow(bp * y, 2) * exp(-pow(bp * y, 2));
    double at = c3 * pow(y, beta) * exp(-alpha * (y - 2));
    double rr = (3.14159265 / 4) * (as + ap + at) * pow(y, 2) / 0.197;
    double dP_dp = rr / 1.01691371;
    if(p>fPCutOff)
      dP_dp = 0;
    assert(dP_dp >= 0);
    // calculate probability density : dProbability/dp
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
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
// Returns the binding energy for a given nucleus.
//____________________________________________________________________________
double EffectiveSF::ReturnBindingEnergy(const Target & target) const
{
  double binding_en;
  if (GetValueFromNuclearMaps(target, fNucRmvE, fRangeNucRmvE, &binding_en) &&
      binding_en > 0) {
    return binding_en;
  }
  return 0;
}
//____________________________________________________________________________
// Returns the fraction of 1p1h events for a given nucleus.  All other events
// are 2p2h.
//____________________________________________________________________________
double EffectiveSF::Returnf1p1h(const Target & target) const
{
  double f1p1h;
  if (GetValueFromNuclearMaps(target, f1p1hMap, fRange1p1hMap, &f1p1h) &&
      f1p1h >= 0 && f1p1h <= 1) {
    return f1p1h;
  }
  return 1;
}
//____________________________________________________________________________
void EffectiveSF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void EffectiveSF::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
// Every parameter for this comes from the config files.
//____________________________________________________________________________
void EffectiveSF::LoadConfig(void)
{
  this->GetParamDef("EjectSecondNucleon2p2h", fEjectSecondNucleon2p2h, false);

  this->GetParamDef("MomentumMax",    fPMax,    1.0);
  this->GetParamDef("MomentumCutOff", fPCutOff, 0.65);
  assert(fPMax > 0 && fPCutOff > 0 && fPCutOff <= fPMax);

  // Find out if Transverse enhancement is enabled to figure out whether to load
  // the 2p2h enhancement parameters.
  this->GetParam("UseElFFTransverseEnhancement", fUseElFFTransEnh );
  if (!fUseElFFTransEnh) {
    LOG("EffectiveSF", pINFO)
        << "Transverse enhancement not used; "
        << "Do not increase the 2p2h cross section.";
  }
  else {
    LoadAllIsotopesForKey("TransEnhf1p1hMod", "EffectiveSF",
                          GetOwnedConfig(), &fTransEnh1p1hMods);
    LoadAllNucARangesForKey("TransEnhf1p1hMod", "EffectiveSF",
                            GetOwnedConfig(), &fRangeTransEnh1p1hMods);
  }

  LoadAllIsotopesForKey("BindingEnergy", "EffectiveSF", GetOwnedConfig(), &fNucRmvE);
  LoadAllNucARangesForKey("BindingEnergy", "EffectiveSF",
                          GetOwnedConfig(), &fRangeNucRmvE);
  LoadAllIsotopesForKey("f1p1h", "EffectiveSF", GetOwnedConfig(), &f1p1hMap);
  LoadAllNucARangesForKey("f1p1h", "EffectiveSF", GetOwnedConfig(), &fRange1p1hMap);

  for (int Z = 1; Z < 140; Z++) {
    for (int A = Z; A < 3*Z; A++) {
      const int pdgc = pdg::IonPdgCode(A, Z);
      double bs, bp, alpha, beta, c1, c2, c3;
      if (GetDoubleKeyPDG("bs", pdgc, GetOwnedConfig(), &bs) &&
          GetDoubleKeyPDG("bp", pdgc, GetOwnedConfig(), &bp) &&
          GetDoubleKeyPDG("alpha", pdgc, GetOwnedConfig(), &alpha) &&
          GetDoubleKeyPDG("beta", pdgc, GetOwnedConfig(), &beta) &&
          GetDoubleKeyPDG("c1", pdgc, GetOwnedConfig(), &c1) &&
          GetDoubleKeyPDG("c2", pdgc, GetOwnedConfig(), &c2) &&
          GetDoubleKeyPDG("c3", pdgc, GetOwnedConfig(), &c3)) {
        vector<double> pars = vector<double>();
        pars.push_back(bs);
        pars.push_back(bp);
        pars.push_back(alpha);
        pars.push_back(beta);
        pars.push_back(c1);
        pars.push_back(c2);
        pars.push_back(c3);
        LOG("EffectiveSF", pINFO)
          << "Nucleus: " << pdgc << " -> using bs =  " << bs << "; bp = "<< bp
          << "; alpha = " << alpha << "; beta = "<<beta<<"; c1 = "<<c1
          <<"; c2 = "<<c2<< "; c3 = " << c3;
        fProbDistParams[pdgc] = pars;
      }
    }
  }
  for(int lowA = 1; lowA < 3 * 140; lowA++) {
    for(int highA = lowA; highA < 3 * 140; highA++) {
      double bs, bp, alpha, beta, c1, c2, c3;
      if (GetDoubleKeyRangeNucA("bs",   lowA, highA, GetOwnedConfig(), &bs)    &&
          GetDoubleKeyRangeNucA("bp",   lowA, highA, GetOwnedConfig(), &bp)    &&
          GetDoubleKeyRangeNucA("alpha",lowA, highA, GetOwnedConfig(), &alpha) &&
          GetDoubleKeyRangeNucA("beta", lowA, highA, GetOwnedConfig(), &beta) &&
          GetDoubleKeyRangeNucA("c1",   lowA, highA, GetOwnedConfig(), &c1) &&
          GetDoubleKeyRangeNucA("c2",   lowA, highA, GetOwnedConfig(), &c2) &&
          GetDoubleKeyRangeNucA("c3",   lowA, highA, GetOwnedConfig(), &c3)) {
        vector<double> pars = vector<double>();
        pars.push_back(bs);
        pars.push_back(bp);
        pars.push_back(alpha);
        pars.push_back(beta);
        pars.push_back(c1);
        pars.push_back(c2);
        pars.push_back(c3);
        LOG("EffectiveSF", pINFO) << "For "<< lowA - 1 <<" < A < " << highA + 1
          <<" -> using bs =  " << bs << "; bp = "<< bp
          << "; alpha = " << alpha << "; beta = "<<beta<<"; c1 = "<<c1
          <<"; c2 = "<<c2<< "; c3 = " << c3;
        fRangeProbDistParams[pair<int, int>(lowA, highA)] = pars;
      }
    }
  }
}
//____________________________________________________________________________
