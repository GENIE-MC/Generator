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

#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "EVGCore/PhysInteractionSelector.h"
#include "EVGCore/EventRecord.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/XSecAlgorithmMap.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"

using std::vector;
using namespace genie;
using namespace genie::units;

//___________________________________________________________________________
PhysInteractionSelector::PhysInteractionSelector() :
InteractionSelectorI("genie::PhysInteractionSelector")
{

}
//___________________________________________________________________________
PhysInteractionSelector::PhysInteractionSelector(string config) :
InteractionSelectorI("genie::PhysInteractionSelector", config)
{

}
//___________________________________________________________________________
PhysInteractionSelector::~PhysInteractionSelector()
{

}
//___________________________________________________________________________
EventRecord * PhysInteractionSelector::SelectInteraction
          (const XSecAlgorithmMap * xscmap, const TLorentzVector & p4) const
{
  if(!xscmap) {
     LOG("InteractionSelector", pERROR)
               << "\n*** NULL XSecAlgorithmMap! Can't select interaction";
     return 0;
  }
  if(xscmap->size() <= 0) {
     LOG("InteractionSelector", pERROR)
              << "\n*** Empty XSecAlgorithmMap! Can't select interaction";
     return 0;
  }

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = 0;
  if (fUseSplines) xssl = XSecSplineList::Instance();

  const InteractionList & ilst = xscmap->GetInteractionList();

  vector<double> xseclist(ilst.size());

  unsigned int i=0;
  InteractionList::const_iterator intliter; // interaction list iter

  for(intliter = ilst.begin(); intliter != ilst.end(); ++intliter) {

     Interaction * interaction = new Interaction(**intliter);
     interaction->GetInitialStatePtr()->SetProbeP4(p4);

     SLOG("InteractionSelector", pDEBUG)
           << "Computing xsec for: \n  " << interaction->AsString();

     // get the cross section for this interaction
     const XSecAlgorithmI * xsec_alg =
                             xscmap->FindXSecAlgorithm(interaction);
     assert(xsec_alg);

     double xsec = 0; // cross section for this interaction

     bool spline_computed = xssl->SplineExists(xsec_alg, interaction);
     bool eval = fUseSplines && spline_computed;
     if (eval) {
           const InitialState & init = interaction->GetInitialState();
           double E = init.GetProbeE(kRfStruckNucAtRest);
           const Spline * spl = xssl->GetSpline(xsec_alg,interaction);
           if(spl->ClosestKnotValueIsZero(E,"-")) xsec = 0;
           else xsec = spl->Evaluate(E);
     } else {
           xsec = xsec_alg->XSec(interaction);
     }
     SLOG("InteractionSelector", pINFO)
       << "\n   " << interaction->AsString() << "\n    --> cross section "
            << (eval ? "[**interpolated**]" : "[**calculated**]") << " = "
               << xsec/cm2 << " cm^2";

     xseclist[i++] = xsec;
     delete interaction;

  } // loop over interaction that can be generated

  // select an interaction

  LOG("InteractionSelector", pINFO)
                        << "Selecting an entry from the Interaction List";
  double xsec_sum  = 0;
  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {
     xsec_sum       += xseclist[iint];
     xseclist[iint]  = xsec_sum;

     SLOG("InteractionSelector", pINFO)
                           << "Sum{xsec}(0->" << iint << ") = " << xsec_sum;
  }
  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->Random1().Rndm();

  LOG("InteractionSelector", pINFO)
               << "Generating Rndm (0. -> max = " << xsec_sum << ") = " << R;

  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {

     SLOG("InteractionSelector", pDEBUG)
                       << "Sum{xsec}(0->" << iint <<") = " << xseclist[iint];

     if( R < xseclist[iint] ) {
       Interaction * selected_interaction = new Interaction (*ilst[iint]);
       selected_interaction->GetInitialStatePtr()->SetProbeP4(p4);

       // set the cross section for the selected interaction (just extract it
       // from the array of summed xsecs rather than recomputing it)
       double xsec_pedestal = (iint > 0) ? xseclist[iint-1] : 0.;
       double xsec = xseclist[iint] - xsec_pedestal;
       assert(xsec>0);

       LOG("InteractionSelector", pNOTICE)
                      << "Selected interaction: \n" << *selected_interaction;

       // bootstrap the event record
       EventRecord * evrec = new EventRecord;
       evrec->AttachInteraction(selected_interaction);
       evrec->SetXSec(xsec);

       return evrec;
     }
  }
  LOG("InteractionSelector", pERROR) << "\nCould not select interaction";
  return 0;
}
//___________________________________________________________________________
void PhysInteractionSelector::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void PhysInteractionSelector::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void PhysInteractionSelector::LoadConfigData(void)
{
  //check whether the user prefers the cross sections to be calculated or
  //evaluated from a spline object constructed at the job initialization
  fUseSplines = fConfig->GetBoolDef("use-stored-xsecs", false);
}
//___________________________________________________________________________
