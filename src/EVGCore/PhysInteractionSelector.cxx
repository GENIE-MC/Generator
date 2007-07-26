//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - January 25, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <vector>
#include <sstream>

#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "EVGCore/PhysInteractionSelector.h"
#include "EVGCore/EventRecord.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionGeneratorMap.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"
#include "Utils/PrintUtils.h"

using std::vector;
using std::ostringstream;
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
     (const InteractionGeneratorMap * igmap, const TLorentzVector & p4) const
{
  if(!igmap) {
     LOG("IntSel", pERROR)
        << "\n*** NULL InteractionGeneratorMap! Can't select interaction";
     return 0;
  }
  if(igmap->size() <= 0) {
     LOG("IntSel", pERROR)
       << "\n*** Empty InteractionGeneratorMap! Can't select interaction";
     return 0;
  }

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = 0;
  if (fUseSplines) xssl = XSecSplineList::Instance();

  const InteractionList & ilst = igmap->GetInteractionList();
  vector<double> xseclist(ilst.size());

  string istate = ilst[0]->InitState().AsString();
  ostringstream msg;
  msg << "Selecting an interaction for the given initial state = "
      << istate << " at E = " << p4.E() << " GeV";

  LOG("IntSel", pNOTICE)
             << utils::print::PrintFramedMesg(msg.str(), 0, '=');
  LOG("IntSel", pNOTICE)
     << "Computing xsecs for all relevant modeled interactions:";

  unsigned int i=0;
  InteractionList::const_iterator intliter = ilst.begin(); 

  for( ; intliter != ilst.end(); ++intliter) {

     Interaction * interaction = new Interaction(**intliter);
     interaction->InitStatePtr()->SetProbeP4(p4);

     SLOG("IntSel", pDEBUG)
           << "Computing xsec for: \n  " << interaction->AsString();

     // get the cross section for this interaction
     const XSecAlgorithmI * xsec_alg =
               igmap->FindGenerator(interaction)->CrossSectionAlg();
     assert(xsec_alg);

     double xsec = 0; // cross section for this interaction

     bool spline_computed = xssl->SplineExists(xsec_alg, interaction);
     bool eval = fUseSplines && spline_computed;
     if (eval) {
           const InitialState & init = interaction->InitState();
           double E = init.ProbeE(kRfHitNucRest);
           const Spline * spl = xssl->GetSpline(xsec_alg,interaction);
           if(spl->ClosestKnotValueIsZero(E,"-")) xsec = 0;
           else xsec = spl->Evaluate(E);
     } else {
           xsec = xsec_alg->Integral(interaction);
     }
     BLOG("IntSel", pNOTICE)
       << interaction->AsString() 
       << " --> xsec " << (eval ? "[**interp**]" : "[**calc**]") 
       << " = " << xsec/cm2 << " cm^2";

     xseclist[i++] = xsec;
     delete interaction;

  } // loop over interaction that can be generated

  // select an interaction

  LOG("IntSel", pINFO)
            << "Selecting an entry from the Interaction List";
  double xsec_sum  = 0;
  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {
     xsec_sum       += xseclist[iint];
     xseclist[iint]  = xsec_sum;

     SLOG("IntSel", pINFO)
             << "Sum{xsec}(0->" << iint << ") = " << xsec_sum;
  }
  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->RndISel().Rndm();

  LOG("IntSel", pINFO)
      << "Generating Rndm (0. -> max = " << xsec_sum << ") = " << R;

  for(unsigned int iint = 0; iint < xseclist.size(); iint++) {

     SLOG("IntSel", pDEBUG)
               << "Sum{xsec}(0->" << iint <<") = " << xseclist[iint];

     if( R < xseclist[iint] ) {
       Interaction * selected_interaction = new Interaction (*ilst[iint]);
       selected_interaction->InitStatePtr()->SetProbeP4(p4);

       // set the cross section for the selected interaction (just extract it
       // from the array of summed xsecs rather than recomputing it)
       double xsec_pedestal = (iint > 0) ? xseclist[iint-1] : 0.;
       double xsec = xseclist[iint] - xsec_pedestal;
       assert(xsec>0);

       LOG("IntSel", pNOTICE)
                << "\nSelected interaction: " << *selected_interaction;

       // bootstrap the event record
       EventRecord * evrec = new EventRecord;
       evrec->AttachSummary(selected_interaction);
       evrec->SetXSec(xsec);

       return evrec;
     }
  }
  LOG("IntSel", pERROR) << "Could not select interaction";
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
  fUseSplines = fConfig->GetBoolDef("UseStoredXSecs", false);
}
//___________________________________________________________________________
