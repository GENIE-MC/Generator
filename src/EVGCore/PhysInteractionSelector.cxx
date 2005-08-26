//____________________________________________________________________________
/*!

\class   genie::PhysInteractionSelector

\brief   Selects interactions to be generated

         Is a concrete implementation of the InteractionSelectorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created January 25, 2005

*/
//____________________________________________________________________________

#include <vector>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "EVGCore/PhysInteractionSelector.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionFilter.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"

using std::vector;
using namespace genie;
using namespace genie::units;

//___________________________________________________________________________
PhysInteractionSelector::PhysInteractionSelector() :
InteractionSelectorI()
{
  fName = "genie::PhysInteractionSelector";

  fEventGeneratorList = 0;
  fInteractionFilter  = 0;
}
//___________________________________________________________________________
PhysInteractionSelector::PhysInteractionSelector(const char * param_set) :
InteractionSelectorI(param_set)
{
  fName = "genie::PhysInteractionSelector";

  this->FindConfig();

  fEventGeneratorList = 0;
  fInteractionFilter  = 0;
}
//___________________________________________________________________________
PhysInteractionSelector::~PhysInteractionSelector()
{

}
//___________________________________________________________________________
void PhysInteractionSelector::SetGeneratorList(
                                          const EventGeneratorList * evglist)
{
  fEventGeneratorList = evglist;
}
//___________________________________________________________________________
void PhysInteractionSelector::SetInteractionFilter(
                                            const InteractionFilter * filter)
{
  fInteractionFilter = filter;
}
//___________________________________________________________________________
Interaction * PhysInteractionSelector::SelectInteraction (
                                       const InitialState & init_state) const
{
  if(!fEventGeneratorList) {
     LOG("InteractionSelector", pERROR)
               << "\n*** NULL Generator List! "
                         << "Can not select interaction for " << init_state;
     return 0;
  }
  if(fEventGeneratorList->size() <= 0) {
     LOG("InteractionSelector", pERROR)
               << "\n*** Empty Generator List! "
                         << "Can not select interaction for " << init_state;
     return 0;
  }

  XSecSplineList * xssl = 0;
  bool used_stored_xsec = this->UseStoredCrossSections();

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  if (used_stored_xsec) xssl = XSecSplineList::Instance();

  vector<double>  xseclist;
  InteractionList intlist;

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  for(evgliter = fEventGeneratorList->begin();
                       evgliter != fEventGeneratorList->end(); ++evgliter) {

     const EventGeneratorI * evgen = *evgliter;

     LOG("InteractionList", pINFO)
            << "Querying [" << evgen->Name() << "/"
                        << evgen->ParamSet() << "] for its InteractionList";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state

     const InteractionListGeneratorI * intlistgen = evgen->IntListGenerator();
     InteractionList * ilst = intlistgen->CreateInteractionList(init_state);

     // cross section algorithm used by this EventGenerator

     const XSecAlgorithmI * xsec_alg = evgen->CrossSectionAlg();

     // loop over all interaction that can be genererated and ask the
     // appropriate cross section algorithm to compute its cross section

     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter) {

         Interaction * interaction = *intliter;

         // get the cross section for this interaction
         SLOG("InteractionList", pINFO)
                              << "Interaction = " << interaction->AsString();

         double xsec = 0; // cross section for this interaction

         bool spline_computed = xssl->SplineExists(xsec_alg, interaction);
         bool eval = used_stored_xsec && spline_computed;

         if (eval) {
           const InitialState & init = interaction->GetInitialState();
           double E = init.GetProbeE(kRfStruckNucAtRest);
           const Spline * spl = xssl->GetSpline(xsec_alg,interaction);
           if(spl->ClosestKnotValueIsZero(E,"-")) xsec = 0;
           else xsec = spl->Evaluate(E);
           SLOG("InteractionList", pINFO)
             << "Cross Section [**interpolated**] = " << xsec/cm2 << " cm^2";
         } else {
           xsec = xsec_alg->XSec(interaction);
           SLOG("InteractionList", pINFO)
               << "Cross Section [**calculated**] = " << xsec/cm2 << " cm^2";
         }

         // save the interaction and its corresponding cross section in the
         // interaction and cross section lists respectivelly.

         intlist.push_back(interaction);
         xseclist.push_back(xsec);

     } // loop over interaction that can be generated
  } // loop over event generators

  // select an interaction

  LOG("InteractionList", pINFO)
                        << "Selecting an entry from the Interaction List";
  double xsec_sum  = 0;
  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {
     xsec_sum       += xseclist[iint];
     xseclist[iint]  = xsec_sum;

     SLOG("InteractionList", pINFO)
                           << "Sum{xsec}(0->" << iint << ") = " << xsec_sum;
  }
  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->Random2().Rndm();

  LOG("InteractionSelector", pINFO)
               << "Generating Rndm (0. -> max = " << xsec_sum << ") = " << R;

  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {

     SLOG("InteractionSelector", pDEBUG)
                       << "Sum{xsec}(0->" << iint <<") = " << xseclist[iint];

     if( R < xseclist[iint] ) {
       Interaction * selected_interaction = new Interaction (*intlist[iint]);

       // set the cross section for the selected interaction (just extract it
       // from the array of summed xsecs rather than recomputing it)
       double xsec_pedestal = (iint > 0) ? xseclist[iint-1] : 0.;
       double xsec = xseclist[iint] - xsec_pedestal;
       assert(xsec>0);
       selected_interaction->SetXSec(xsec);

       LOG("InteractionSelector", pINFO)
                      << "Selected interaction: \n" << *selected_interaction;
       return selected_interaction;
     }
  }
  LOG("InteractionSelector", pERROR)
                      << "\nCould not select interaction for " << init_state;
  return 0;
}
//___________________________________________________________________________
bool PhysInteractionSelector::UseStoredCrossSections(void) const
{
//check whether the user prefers the cross sections to be calculated or
//evaluated from a spline object constructed at the job initialization

  return fConfig->Exists("use-stored-xsecs") ?
                                fConfig->GetBool("use-stored-xsecs") : false;
}
//___________________________________________________________________________
