//____________________________________________________________________________
/*!

\class    genie::NuclMomentumGenerator

\brief    Holds a list of all nucleon momentum probability distributions that
          have been generated at a MC job. Can pick up the one requested (or
          build it if it doesn't already exists) and use it to generate nucleon
          momenta.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#include <sstream>
#include <iostream>

#include <TMath.h>
#include <TH1D.h>
#include <TVector3.h>

#include "Conventions/Constants.h"
#include "Interaction/Target.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclMomentumModelI.h"
#include "Nuclear/NuclMomentumGenerator.h"
#include "Numerical/RandomGen.h"

using std::ostringstream;
using std::cout;
using std::endl;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuclMomentumGenerator * NuclMomentumGenerator::fInstance = 0;
//____________________________________________________________________________
NuclMomentumGenerator::NuclMomentumGenerator()
{
  fInstance = 0;
  fCurrProbDistribution = 0;
}
//____________________________________________________________________________
NuclMomentumGenerator::~NuclMomentumGenerator()
{
  cout << "NuclMomentumGenerator singleton dtor: "
                  << "Deleting all nucleon momentum distributions" << endl;

  map<string, TH1D*>::iterator it;
  for(it=fProbDistributionMap.begin(); it!=fProbDistributionMap.end(); ++it) {
    TH1D * histo = it->second;
    //if(histo) delete histo; ?
    histo=0;
  }
  fProbDistributionMap.clear();
  fInstance = 0;
  fCurrProbDistribution = 0;
}
//____________________________________________________________________________
NuclMomentumGenerator * NuclMomentumGenerator::Instance()
{
  if(fInstance == 0) {
    LOG("NucleonP", pINFO) << "NuclMomentumGenerator late initialization";

    static NuclMomentumGenerator::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new NuclMomentumGenerator;
  }
  return fInstance;
}
//____________________________________________________________________________
void NuclMomentumGenerator::UseProbDistribution(
                        const NuclMomentumModelI * model, const Target & tgt)
{
  assert(model);
  assert(tgt.StruckNucleonIsSet());

  fCurrProbDistribution = 0;
  string key = this->BuildProbDistributionKey(model, tgt);

  map<string, TH1D*>::iterator it = fProbDistributionMap.find(key);

  if(it != fProbDistributionMap.end()) {
    fCurrProbDistribution = it->second;
  } else {
    LOG("NucleonP", pNOTICE)
      << "Requested prob=f(Pnucl) distribution not built yet. Building now.";
    LOG("NucleonP", pNOTICE)
      << "Running alg: "<< model->Id().Key() << " on tgt: "<< tgt.AsString();

    TH1D * prob_distribution = model->ProbabilityDistribution(tgt);
    assert(prob_distribution);
    fProbDistributionMap.insert(
               map<string, TH1D*>::value_type(key,prob_distribution));
    fCurrProbDistribution = prob_distribution;
  }
}
//____________________________________________________________________________
double NuclMomentumGenerator::Probability(const TVector3 & p) const
{
  return this->Probability(p.Mag());
}
//____________________________________________________________________________
double NuclMomentumGenerator::Probability(double p) const
{
  if(!fCurrProbDistribution) {
    LOG("NucleonP", pFATAL)
         << "Null nucleon momentum probability distribution";
  }
  assert(fCurrProbDistribution);

  int bin = this->fCurrProbDistribution->FindBin(p);
  double prob = (double) this->fCurrProbDistribution->GetBinContent(bin);
  LOG("NucleonP", pINFO) << "Prob(|p,nucleon| = " << p << ") = " << prob;
  return prob;
}
//____________________________________________________________________________
double NuclMomentumGenerator::RandomMomentum(void) const
{
  if(!fCurrProbDistribution) {
    LOG("NucleonP", pFATAL)
         << "Null nucleon momentum probability distribution";
  }
  assert(fCurrProbDistribution);

  double p = fCurrProbDistribution->GetRandom();
  LOG("NucleonP", pINFO) << "|p,nucleon| = " << p;
  return p;
}
//____________________________________________________________________________
TVector3 NuclMomentumGenerator::RandomMomentum3(void) const
{
  double p = this->RandomMomentum();

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  TVector3 p3(px,py,pz);

  return p3;
}
//____________________________________________________________________________
string NuclMomentumGenerator::BuildProbDistributionKey(
                         const NuclMomentumModelI * model, const Target & tgt)
{
  ostringstream key;
  key << model->Id().Key() << "/" << tgt.AsString();
  return key.str();
}
//____________________________________________________________________________

