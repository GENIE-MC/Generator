//____________________________________________________________________________
/*!

\class    genie::NuclearPDistribution

\brief    Describes a nucleon momentum probability distribution (constructed
          from the attached NuclearPDistributionModelI) and can act as a
          nucleon momentum generator.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearPDistribution.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuclearPDistribution::NuclearPDistribution()
{
  this->Init();
}
//____________________________________________________________________________
NuclearPDistribution::NuclearPDistribution(const NuclearPDistribution & npd)
{
  if (this->fProbDistribution) delete this->fProbDistribution;

  this->Init();

  this->fProbModel = npd.fProbModel;
  this->fProbDistribution = new TH1D( *(npd.fProbDistribution) );
}
//____________________________________________________________________________
NuclearPDistribution::~NuclearPDistribution()
{
  if (this->fProbDistribution) delete this->fProbDistribution;
}
//____________________________________________________________________________
void NuclearPDistribution::AttachModel(const NuclearPDistributionModelI * model)
{
  LOG("NucleonP", pDEBUG) << "Initializing & attaching model";

  if (this->fProbDistribution) delete this->fProbDistribution;
  this->fProbDistribution = 0;

  this->fProbModel = model;
}
//____________________________________________________________________________
void NuclearPDistribution::BuildProbDistribution(const Target & target)
{
  LOG("NucleonP", pDEBUG) << "Building probability distributions";

  this->fProbDistribution = fProbModel->ProbabilityDistribution(target);
}
//____________________________________________________________________________
double NuclearPDistribution::Probability(const TVector3 & p) const
{
  return this->Probability( p.Mag() );
}
//____________________________________________________________________________
double NuclearPDistribution::Probability(double p) const
{
  if (this->fProbDistribution) {
     int bin = this->fProbDistribution->FindBin(p);
     return (double) this->fProbDistribution->GetBinContent(bin);
  } else {
     LOG("NucleonP", pERROR) << "Probability distribution is not built yet";
  }
  return -1;
}
//____________________________________________________________________________
double NuclearPDistribution::RandomMomentum(void) const
{
  if (this->fProbDistribution) {
     double p = this->fProbDistribution->GetRandom();
     LOG("NucleonP", pINFO) << "|p,nucleon| = " << p;
     return p;
  } else {
     LOG("NucleonP", pERROR) << "Probability distribution is not built yet";
  }
  return 0;
}
//____________________________________________________________________________
TVector3 NuclearPDistribution::RandomMomentum3(void) const
{
  double p = this->RandomMomentum();

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->Random2().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->Random2().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  TVector3 p3(px,py,pz);

  return p3;
}
//____________________________________________________________________________
void NuclearPDistribution::Init(void)
{
  this->fProbModel        = 0;
  this->fProbDistribution = 0;
}
//____________________________________________________________________________



