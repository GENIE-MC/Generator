//____________________________________________________________________________
/*!

\class    genie::NuclearPDistribution

\brief    Describes a Nucleon Momentum Probability Distribution.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearPDistribution.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuclearPDistribution::NuclearPDistribution()
{
  Init();
}
//____________________________________________________________________________
NuclearPDistribution::NuclearPDistribution(const NuclearPDistribution & npd)
{
  Init();

  fProbModel = npd.fProbModel;

  fProbDistribution = new TH1D( *(npd.fProbDistribution) );
}
//____________________________________________________________________________
NuclearPDistribution::~NuclearPDistribution()
{
  if (fProbDistribution) delete fProbDistribution;
}
//____________________________________________________________________________
void NuclearPDistribution::AttachModel(const NuclearPDistributionModelI * model)
{
  LOG("NucleonP", pDEBUG) << "Initializing & attaching model";

  if (fProbDistribution) delete fProbDistribution;

  fProbDistribution = 0;

  fProbModel    = model;
}
//____________________________________________________________________________
void NuclearPDistribution::BuildProbDistribution(const Target & target)
{
  LOG("NucleonP", pDEBUG) << "Building probability distributions";

  fProbDistribution = fProbModel->ProbabilityDistribution(target);
}
//____________________________________________________________________________
double NuclearPDistribution::Probability(const TVector3 & p) const
{
  return Probability( p.Mag() );
}
//____________________________________________________________________________
double NuclearPDistribution::Probability(double p) const
{
  if (fProbDistribution) {

     int bin = fProbDistribution->FindBin(p);
     return (double) fProbDistribution->GetBinContent(bin);

  } else {
     LOG("NucleonP", pERROR) << "Probability distribution is not built yet";
  }

  return -1;
}
//____________________________________________________________________________
double NuclearPDistribution::RandomMomentum(void) const
{
  if (fProbDistribution) {

     double p = fProbDistribution->GetRandom();

     SLOG("NucleonP", pINFO) << "|p,nucleon| = " << p;
     
     return p;

  } else {
     LOG("NucleonP", pERROR) << "Probability distribution is not built yet";
  }

  return 0;
}
//____________________________________________________________________________
TVector3 NuclearPDistribution::RandomMomentum3(void) const
{
  TVector3 p(0, 0, RandomMomentum()); // fix in z direction

  // get Euler angles for a random rotation of the momentum 3-vector
  
  RandomGen * rnd = RandomGen::Instance();

  double theta = 2*kPi * rnd->Random2().Rndm();
  double fi    =   kPi * rnd->Random2().Rndm();
  double psi   = 2*kPi * rnd->Random2().Rndm();

  SLOG("NucleonP", pINFO) << "Euler angles: " <<
                 "theta = " << theta << ", fi = " << fi << ", psi = " << psi;
   
  TRotation random_rotation;

  random_rotation.SetToIdentity();
  random_rotation.SetXEulerAngles(theta, fi, psi);

  p.Transform(random_rotation);
  
  return p;
}
//____________________________________________________________________________
void NuclearPDistribution::Init(void)
{
  fProbModel        = 0;
  fProbDistribution = 0;
}
//____________________________________________________________________________



