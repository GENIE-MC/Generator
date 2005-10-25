//____________________________________________________________________________
/*!

\class   genie::WeightCalculator

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#include "Algorithm/AlgCmp.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EventRecord.h"
#include "ReWeight/WeightCalculator.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
WeightCalculator::WeightCalculator()
{

}
//___________________________________________________________________________
WeightCalculator::~WeightCalculator()
{

}
//___________________________________________________________________________
void WeightCalculator::OldCrossSectionModel(const MCModel & model)
{
  fOldMCModel.Copy(model);
}
//___________________________________________________________________________
void WeightCalculator::NewCrossSectionModel(const MCModel & model)
{
  fNewMCModel.Copy(model);
}
//___________________________________________________________________________
double WeightCalculator::ReWeight(const EventRecord & event)
{
  // get event summary (Interaction)
  Interaction * interaction = event.GetInteraction();

  LOG("ReWeight", pDEBUG) << "Computing new weight for: \n" << *interaction;

  // get the appropriate xsec alhorithms
  const XSecAlgorithmI * old_alg = fOldMCModel.XSecAlg(interaction);
  const XSecAlgorithmI * new_alg = fNewMCModel.XSecAlg(interaction);

  // if one of the models is null then warn and set the weight to 0
  if(!old_alg || !new_alg) {
     LOG("ReWeight", pWARN)
         << "Got a NULL XSec Model for: \n" << interaction->AsString();
     return 0.;
  }

  // if the old and new xsec models are identical, then do not
  // bother computing a weight
  if( old_alg->Compare(new_alg) == kAlgCmpIdentical ) {

    LOG("ReWeight", pDEBUG) 
        << "Same old/new xsec models for the given process. Weight = 1.";
    return 1.;
  }

  double old_xsec   = old_alg->XSec(interaction);
  double new_xsec   = new_alg->XSec(interaction);
  double old_weight = event.GetWeight();
  double new_weight = old_weight * (new_xsec/old_xsec);

  LOG("ReWeight", pINFO)
              << "Event weight: " << old_weight << " ---> " << new_weight;

  return new_weight;
}
//___________________________________________________________________________

