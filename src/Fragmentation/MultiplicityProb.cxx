//____________________________________________________________________________
/*!

\class    genie::MultiplicityProb

\brief    Describes a Multiplicity Probability Distribution.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#include "Fragmentation/MultiplicityProb.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
MultiplicityProb::MultiplicityProb()
{
  Init();
}
//____________________________________________________________________________
MultiplicityProb::MultiplicityProb(const MultiplicityProb & mpd)
{
  Init();

  fMultProbModel    = mpd.fMultProbModel;

  fProbDistribution = new TH1D( *(mpd.fProbDistribution) );
}
//____________________________________________________________________________
MultiplicityProb::~MultiplicityProb()
{
  if (fProbDistribution) delete fProbDistribution;
}
//____________________________________________________________________________
void MultiplicityProb::AttachModel(const MultiplicityProbModelI * model)
{
  LOG("MultProb", pINFO) << "Initializing & attaching model";

  if (fProbDistribution) delete fProbDistribution;

  fProbDistribution = 0;
  
  fMultProbModel    = model;
}
//____________________________________________________________________________
void MultiplicityProb::BuildProbDistribution(const Interaction * interaction)
{
  LOG("MultProb", pINFO) << "Building probability distributions";

  fProbDistribution = fMultProbModel->ProbabilityDistribution(interaction);
}
//____________________________________________________________________________
double MultiplicityProb::Probability(int n) const
{
  if (fProbDistribution) {

     int bin = fProbDistribution->FindBin(n);
     return (double) fProbDistribution->GetBinContent(bin);

  } else {
     LOG("MultProb", pERROR) << "Probability distribution is not built yet";
  }

  return -1;
}
//____________________________________________________________________________
unsigned int MultiplicityProb::RandomMultiplicity(
                                     unsigned int min, unsigned int max) const
{
  assert(min < max);
  
  if (fProbDistribution) {

     unsigned int mult = (unsigned int) fProbDistribution->GetRandom();

     // re-try if multiplicity is not at the requested range     
     if(mult < min || mult > max) return RandomMultiplicity(min, max);

     return mult;
     
  } else {
     LOG("MultProb", pERROR) << "Probability distribution is not built yet";
  }

  return 0;
}
//____________________________________________________________________________
void MultiplicityProb::Init(void) 
{
  fMultProbModel    = 0;

  fProbDistribution = 0;
}
//____________________________________________________________________________



