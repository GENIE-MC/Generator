//____________________________________________________________________________
/*!

\class   genie::GEVGDriver

\brief   Minimal interface object for generating neutrino interactions for
         a given initial state.

         When GMC is used, a GEVGDriver object list is assembled for all possible
         initial states (corresponding to combinations of all neutrino types
         -declared by the input GFluxI- and all target nuclei types -found
         in the input ROOT geometry-.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 06, 2004

*/
//____________________________________________________________________________

#include <cassert>

#include <TSystem.h>
#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGCore/EventRecord.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EGResponsibilityChain.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/ToyInteractionSelector.h"
#include "EVGCore/PhysInteractionSelector.h"
#include "EVGCore/EventGeneratorListAssembler.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "EVGCore/InteractionFilter.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::controls;

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const GEVGDriver & driver)
 {
   driver.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
GEVGDriver::GEVGDriver()
{
  this->Initialize();
  this->Configure();
}
//___________________________________________________________________________
GEVGDriver::~GEVGDriver()
{
  if (fNuclTarget)  delete fNuclTarget;
  if (fFilter)      delete fFilter;
  if (fIntSelector) delete fIntSelector;
}
//___________________________________________________________________________
void GEVGDriver::SetInitialState(const InitialState & init_state)
{
  int nu_pdgc = init_state.GetProbePDGCode();
  int Z       = init_state.GetTarget().Z();
  int A       = init_state.GetTarget().A();

  this->SetInitialState(nu_pdgc, Z, A);
}
//___________________________________________________________________________
void GEVGDriver::SetInitialState(int nu_pdgc, int Z, int A)
{
  if (fNuclTarget) delete fNuclTarget;

  fNuPDG      = nu_pdgc;
  fNuclTarget = new Target(Z, A);

  assert( fNuclTarget->IsValidNucleus() || fNuclTarget->IsFreeNucleon() );

  int tgtpdg = fNuclTarget->PDGCode();
  LOG("GEVGDriver", pINFO) << "Set neutrino PDG-code:......... " << fNuPDG;
  LOG("GEVGDriver", pINFO) << "Set nuclear target PDG-code::.. " << tgtpdg;
}
//___________________________________________________________________________
void GEVGDriver::SetFilter(const InteractionFilter & filter)
{
// Sets an InteractionFilter that can supress entries of the InteractionList
// from being selected. This is to be used when one is interested in some
// event classes only (eg QEL CC) and wants to suppress the generation of
// other event types without having to touch the XML configuration files.

  if (fFilter) delete fFilter;

  fFilter = new InteractionFilter(filter);

  if (fFilter && fIntSelector) fIntSelector->SetInteractionFilter(fFilter);
}
//___________________________________________________________________________
void GEVGDriver::FilterUnphysical(bool on_off)
{
  LOG("GEVGDriver", pNOTICE)
        << "*** Filtering unphysical events is turned "
                     << utils::print::BoolAsIOString(on_off) << " ***\n";
  fFilterUnphysical = on_off;
}
//___________________________________________________________________________
void GEVGDriver::Initialize(void)
{
  fCurrentRecord = 0;
  fNuclTarget    = 0;
  fEvGenList     = 0;
  fIntSelector   = 0;
  fChain         = 0;
  fFilter        = 0;
  fUseSplines    = false;
  fNRecLevel     = 0;

  // Default driver behaviour is to filter out unphysical events,
  // Set this to false to get them if needed, but be warned that the event
  // record for unphysical events might be incomplete depending on the
  // processing step that event generation was stopped.
  fFilterUnphysical = true;
}
//___________________________________________________________________________
void GEVGDriver::Configure(void)
{
  LOG("GEVGDriver", pNOTICE) << "Configuring a GEVGDriver object";

  // figure out which list of event generators to use from the $GEVGL
  // environmental variable (use "Default") if the variable is not set.

  LOG("GEVGDriver", pNOTICE) << "Creating the `Event Generator List`";

  string evgl = (gSystem->Getenv("GEVGL") ?
                               gSystem->Getenv("GEVGL") : "Default");
  LOG("GEVGDriver", pNOTICE) << "Specified Event Generator List = " << evgl;

  EventGeneratorListAssembler evglist_assembler(evgl.c_str());
  fEvGenList = evglist_assembler.AssembleGeneratorList();

  LOG("GEVGDriver", pNOTICE) << "Creating the `Generator Chain of Responsibility`";

  if(fChain) delete fChain;
  fChain = new EGResponsibilityChain;
  fChain->SetGeneratorList(fEvGenList);

  LOG("GEVGDriver", pNOTICE) << "Creating an `Interaction Selector`";

  if(fIntSelector) delete fIntSelector;
  fIntSelector = new PhysInteractionSelector("Default");
  fIntSelector->SetGeneratorList(fEvGenList);
}
//___________________________________________________________________________
EventRecord * GEVGDriver::GenerateEvent(const TLorentzVector & nu4p)
{
  //-- Build initial state information from inputs
  LOG("GEVGDriver", pINFO) << "Creating an init state from the user inputs";
  InitialState init_state(*fNuclTarget, fNuPDG);
  init_state.SetProbeP4(nu4p);

  //-- Select the interaction to be generated (amongst the entries of the
  //   InteractionList assembled by the EventGenerators) and bootstrap the
  //   event record
  LOG("GEVGDriver", pINFO)
               << "Selecting an Interaction & Bootstraping the EventRecord";
  fCurrentRecord = fIntSelector->SelectInteraction(init_state);

  assert(fCurrentRecord); // abort if no interaction could be selected!

  //-- Get a ptr to the interaction summary
  LOG("GEVGDriver", pDEBUG)
                     << "Getting the interaction from the new event-record";
  Interaction * interaction = fCurrentRecord->GetInteraction();

  //-- Find the appropriate concrete EventGeneratorI implementation
  //   for generating this event.
  //
  //   The right EventGeneratorI will be selecting by iterating over the
  //   entries of the EventGeneratorList and compare the interaction
  //   against the ValidityContext declared by each EventGeneratorI
  //
  //   (note: use of the 'Chain of Responsibility' Design Pattern)
  LOG("GEVGDriver", pINFO) << "Finding an appropriate EventGenerator";
  const EventGeneratorI * evgen = fChain->FindGenerator(interaction);

  assert(evgen);

  //-- Generate the selected event
  //
  //   The selected EventGeneratorI subclass will start processing the
  //   event record (by sequentially asking each entry in its list of
  //   EventRecordVisitorI subclasses to visit and process the record).
  //   Most of the actual event generation takes place in this step.
  //
  //   (note: use of the 'Visitor' Design Pattern)
  LOG("GEVGDriver", pINFO) << "Generating Event:";
  evgen->ProcessEventRecord(fCurrentRecord);

  //-- If the user requested that unphysical events should be returned too,
  //   return the event record here

  if (!fFilterUnphysical) return fCurrentRecord;

  //-- Check whether the generated event is unphysical (eg Pauli-Blocked,
  //   etc). If the unphysical it will not be returned. This method would
  //   clean-up the failed event and it will enter into a recursive mode,
  //   calling itself so as to regenerate the event. It will exit the
  //   recursive mode when a physical event is generated
  bool unphys = fCurrentRecord->IsUnphysical();
  if(unphys) {
     LOG("GEVGDriver", pWARN) << "I generated an unphysical event!";
     delete fCurrentRecord;
     fCurrentRecord = 0;
     fNRecLevel++; // increase the nested level counter

     if(fNRecLevel<=kRecursiveModeMaxDepth) {
         LOG("GEVGDriver", pWARN) << "Attempting to regenerate the event.";
         return this->GenerateEvent(nu4p);
     } else {
        LOG("GEVGDriver", pFATAL)
             << "Could not produce a physical event after "
                      << kRecursiveModeMaxDepth << " attempts - Aborting!";
        assert(false);
     }
  }
  fNRecLevel = 0;

  //-- Return a successfully generated event to the user
  return fCurrentRecord; // The client 'adopts' the event record
}
//___________________________________________________________________________
double GEVGDriver::SumCrossSection(const TLorentzVector & nup4)
{
  LOG("GEVGDriver", pDEBUG) << "Computing the cross section sum";

  double xsec_sum = 0;

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = XSecSplineList::Instance();

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  // build the initial state
  if(!this->IsValidInitState()) {
    LOG("GEVGDriver", pERROR) << "Invalid initial state"; return 0;
  }
  InitialState init_state(*fNuclTarget, fNuPDG);

  // loop over all EventGenerator objects used in the current job

  for(evgliter = fEvGenList->begin();
                             evgliter != fEvGenList->end(); ++evgliter) {
     const EventGeneratorI * evgen = *evgliter; // current EventGenerator

     LOG("GEVGDriver", pNOTICE)
            << "Querying [" << evgen->Id().Key()
                                           << "] for its InteractionList";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state
     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(init_state);

     //no point to go on if the list is NULL - continue to next iteration
     if(!ilst) continue;

     // cross section algorithm used by this EventGenerator
     const XSecAlgorithmI * xsec_alg = evgen->CrossSectionAlg();

     // loop over all interaction that can be genererated and ask the
     // appropriate cross section algorithm to compute its cross section

     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter) {
         Interaction * interaction = *intliter;

         // set the input 4-momentum
         interaction->GetInitialStatePtr()->SetProbeP4(nup4);

         // get the cross section for this interaction
         SLOG("GEVGDriver", pNOTICE)
                          << "Interaction = " << interaction->AsString();

         double xsec = 0;
         bool spline_exists = xssl->SplineExists(xsec_alg, interaction);
         if (spline_exists && fUseSplines) {
             double E = interaction->GetInitialState().GetProbeE(kRfStruckNucAtRest);
             xsec = xssl->GetSpline(xsec_alg,interaction)->Evaluate(E);
         } else
             xsec = xsec_alg->XSec(interaction);

         xsec_sum += xsec;
         SLOG("GEVGDriver", pNOTICE)
                  << "Cross Section = " << (xsec/units::cm2) << " cm2";
     } // loop over interaction that can be generated by this generator

     delete ilst;
     ilst = 0;
  } // loop over event generators

  return xsec_sum;
}
//___________________________________________________________________________
double GEVGDriver::MaxCrossSection(const TLorentzVector & nup4)
{
  LOG("GEVGDriver", pDEBUG)
          << "Computing the cross section for the max cross section process";

  double max_xsec = 0;

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = XSecSplineList::Instance();

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  // build the initial state
  if(!this->IsValidInitState()) {
    LOG("GEVGDriver", pERROR) << "Invalid initial state"; return 0;
  }
  InitialState init_state(*fNuclTarget, fNuPDG);

  // loop over all EventGenerator objects used in the current job

  for(evgliter = fEvGenList->begin();
                             evgliter != fEvGenList->end(); ++evgliter) {
     const EventGeneratorI * evgen = *evgliter; // current EventGenerator

     LOG("GEVGDriver", pINFO)
            << "Querying [" << evgen->Id().Key()
                                          << "] for its InteractionList";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state
     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(init_state);

     //no point to go on if the list is NULL - continue to next iteration
     if(!ilst) continue;

     // cross section algorithm used by this EventGenerator
     const XSecAlgorithmI * xsec_alg = evgen->CrossSectionAlg();

     // loop over all interaction that can be genererated and ask the
     // appropriate cross section algorithm to compute its cross section

     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter) {
         Interaction * interaction = *intliter;

         // set the input 4-momentum
         interaction->GetInitialStatePtr()->SetProbeP4(nup4);

         // get the cross section for this interaction
         SLOG("GEVGDriver", pNOTICE)
                          << "Interaction = " << interaction->AsString();

         double xsec = 0;
         bool spline_exists = xssl->SplineExists(xsec_alg, interaction);
         if (spline_exists && fUseSplines) {
             double E = interaction->GetInitialState().GetProbeE(kRfStruckNucAtRest);
             xsec = xssl->GetSpline(xsec_alg,interaction)->Evaluate(E);
         } else
             xsec = xsec_alg->XSec(interaction);

         max_xsec = TMath::Max(max_xsec, xsec);
         SLOG("GEVGDriver", pNOTICE)
                  << "Cross Section = " << (xsec/units::cm2) << " cm2";
     } // loop over interaction that can be generated by this generator

     delete ilst;
     ilst = 0;
  } // loop over event generators

  return max_xsec;
}
//___________________________________________________________________________
void GEVGDriver::UseSplines(void)
{
  // build the initial state
  if(!this->IsValidInitState()) {
    LOG("GEVGDriver", pERROR) << "Invalid initial state"; return;
  }
  InitialState init_state(*fNuclTarget, fNuPDG);

  fUseSplines = true;
  XSecSplineList * xsl = XSecSplineList::Instance();

  // if the user wants to use splines, make sure that all the splines needed
  // have been computed or loaded
  if(fUseSplines) {
     EventGeneratorList::const_iterator evgliter;
     InteractionList::iterator          iliter;

     // loop over all EventGenerator objects used in the current job
     for(evgliter = fEvGenList->begin();
                                evgliter != fEvGenList->end(); ++evgliter) {

       const EventGeneratorI * evgen = *evgliter;

       // ask the event generator to produce a list of all interaction it can
       // generate for the input initial state
       const InteractionListGeneratorI * ilgen = evgen->IntListGenerator();
       InteractionList * ilst = ilgen->CreateInteractionList(init_state);
       if(!ilst) continue;

       // total cross section algorithm used by the current EventGenerator
       const XSecAlgorithmI * alg = evgen->CrossSectionAlg();

       // loop over all interaction that can be genererated
       for(iliter = ilst->begin(); iliter != ilst->end(); ++iliter) {

           Interaction * interaction = *iliter;
           bool spl_exists = xsl->SplineExists(alg, interaction);

           fUseSplines = fUseSplines && spl_exists;

           if(!spl_exists) {
              LOG("GEVGDriver", pWARN)
                    << "\nAt least a spline does not exist - "
                                     << "Reverting back to not using splines";
              return;
           }
       } // loop over interaction that can be generated by this generator

     delete ilst;
     ilst = 0;

     } // loop over event generators
  }//use-splines?
}
//___________________________________________________________________________
void GEVGDriver::CreateSplines(bool useLogE)
{
  LOG("GEVGDriver", pINFO)
       << "\nCreating Cross Section Splines with UseLogE = "
                                               << ((useLogE) ? "ON" : "OFF");
  // build the initial state
  if(!this->IsValidInitState()) {
    LOG("GEVGDriver", pERROR) << "Invalid initial state"; return;
  }
  InitialState init_state(*fNuclTarget, fNuPDG);

  XSecSplineList * xsl = XSecSplineList::Instance();
  xsl->SetLogE(useLogE);

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  // loop over all EventGenerator objects used in the current job

  for(evgliter = fEvGenList->begin();
                               evgliter != fEvGenList->end(); ++evgliter) {

     const EventGeneratorI * evgen = *evgliter;

     LOG("GEVGDriver", pNOTICE)
            << "Querying [ " << evgen->Id().Key()
                                            << "] for its InteractionList";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state
     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(init_state);
     if(!ilst) continue;

     // total cross section algorithm used by the current EventGenerator
     const XSecAlgorithmI * alg = evgen->CrossSectionAlg();

     // get the energy range of the spline from the EventGenerator
     // validity context
     double Emin = TMath::Max(0.01,evgen->ValidityContext().Emin());
     double Emax = evgen->ValidityContext().Emax();

     // loop over all interaction that can be genererated and ask the
     // appropriate cross section algorithm to compute its cross section

     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter) {

         Interaction * interaction = *intliter;

         // create a cross section spline for this interaction & store
         LOG("GEVGDriver", pNOTICE)
                    << "\nCreating xsec spline for \n" << *interaction;

         // only create the spline if it does not already exists
         bool spl_exists = xsl->SplineExists(alg, interaction);
         if(!spl_exists) {
             LOG("GEVGDriver", pINFO) << "Computing spline knots";
             xsl->CreateSpline(alg, interaction, 40, Emin, Emax);
         } else {
             LOG("GEVGDriver", pNOTICE) << "Spline is already loaded - skipping";
         }

     } // loop over interaction that can be generated by this generator

     delete ilst;
     ilst = 0;

  } // loop over event generators

  LOG("GEVGDriver", pINFO) << *xsl; // print list of splines

  fUseSplines = true;
}
//___________________________________________________________________________
bool GEVGDriver::IsValidInitState(void) const
{
  if(!fNuclTarget) return false;

  bool isnu = pdg::IsNeutrino(fNuPDG) || pdg::IsAntiNeutrino(fNuPDG);
  if(!isnu) return false;

  return true;
}
//___________________________________________________________________________
void GEVGDriver::Print(ostream & stream) const
{
  stream
    << "\n\n *********************** GEVGDriver ***************************";

  if( this->IsValidInitState() ) {
    int tgtpdg = fNuclTarget->PDGCode();
    stream << "\n  |---o Neutrino PDG-code .........: " << fNuPDG;
    stream << "\n  |---o Nuclear Target PDG-code ...: " << tgtpdg;
  } else
    stream << "\n  |---o *** The initial state wasn't defined properly ***";

  if (fFilter) {
    stream << "\n  |---o An InteractionFilter is being used: " << *fFilter;
  }
  stream << "\n  |---o Using cross section splines is turned "
                                << utils::print::BoolAsIOString(fUseSplines);
  stream << "\n  |---o Filtering unphysical events is turned "
                          << utils::print::BoolAsIOString(fFilterUnphysical);
  stream << "\n *********************************************************\n";
}
//___________________________________________________________________________
