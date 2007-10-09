//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - August 06, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <sstream>

#include <TSystem.h>
#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGCore/EventRecord.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/ToyInteractionSelector.h"
#include "EVGCore/PhysInteractionSelector.h"
#include "EVGCore/EventGeneratorListAssembler.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "EVGCore/InteractionGeneratorMap.h"
#include "EVGCore/RunningThreadInfo.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/XSecSplineList.h"
#include "Utils/PrintUtils.h"

using std::ostringstream;

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
  this->Init();
}
//___________________________________________________________________________
GEVGDriver::~GEVGDriver()
{
  this->CleanUp();
}
//___________________________________________________________________________
void GEVGDriver::FilterUnphysical(const TBits & unphysmask)
{
  *fFiltUnphysMask = unphysmask;

  LOG("GEVGDriver", pNOTICE) 
    << "Set unphysical event mask (" << GHepFlags::NFlags() << "->0) = " 
    << *fFiltUnphysMask;
}
//___________________________________________________________________________
void GEVGDriver::Init(void)
{
  // initial state for which this driver is configured
  fInitState = 0;

  // current event record (ownership is transfered at GenerateEvent())
  fCurrentRecord = 0;

  // list of Event Generator objects loaded into the driver
  fEvGenList = 0;

  // interaction selector
  fIntSelector = 0;

  // flag instructing the driver whether, for each interaction, to compute
  // cross section by running their corresponding XSecAlgorithm or by 
  // evaluating their corresponding xsec spline
  fUseSplines = false;

  // 'depth' counter when entering a recursive mode to re-generate a failed/
  // unphysical event - the driver is not allowed to go into arbitrarily large
  // depths
  fNRecLevel = 0;

  // an "interaction" -> "generator" associative contained built for all
  // simulated interactions (from the loaded Event Generators and for the 
  // input initial state)
  fIntGenMap = 0;

  // A spline describing the sum of all interaction cross sections given an
  // initial state (the init state with which this driver was configured).
  // Create it using the CreateXSecSumSpline() method
  // The sum of all interaction cross sections is used, for example, by
  // GMCJDriver for selecting an initial state.
  fXSecSumSpl = 0;

  // Default driver behaviour is to filter out unphysical events
  // If needed, set the fFiltUnphysMask bitfield to get pre-selected types of 
  // unphysical events (just set to 1 the bit you want ignored from the check).
  // You can do that either by calling the FilterUnphysical() or setting the 
  // GUNPHYSMASK env.var) but be warned that the event record for unphysical 
  // events might be incomplete depending on the processing step that event 
  // generation was stopped.

  fFiltUnphysMask = new TBits(GHepFlags::NFlags());
  fFiltUnphysMask->ResetAllBits(false); 

  if (gSystem->Getenv("GUNPHYSMASK")) {
     unsigned int i=0;
     const char * bitfield = gSystem->Getenv("GUNPHYSMASK");

     while(bitfield[i]) {
        bool flag = (bitfield[i]=='1');

        if(i<GHepFlags::NFlags()) fFiltUnphysMask->SetBitNumber(i,flag);
        i++;
     }
  }
  LOG("GEVGDriver", pNOTICE) 
    << "Initializing unphysical event mask (" << GHepFlags::NFlags() 
    << "->0) = " << *fFiltUnphysMask;
}
//___________________________________________________________________________
void GEVGDriver::CleanUp(void)
{
  if (fInitState)        delete fInitState;
  if (fEvGenList)        delete fEvGenList;
  if (fIntSelector)      delete fIntSelector;
  if (fIntGenMap)        delete fIntGenMap;
  if (fXSecSumSpl)       delete fXSecSumSpl;
  if (fFiltUnphysMask)   delete fFiltUnphysMask;
}
//___________________________________________________________________________
void GEVGDriver::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void GEVGDriver::Configure(int nu_pdgc, int Z, int A)
{
  Target target(Z, A);
  InitialState init_state(target, nu_pdgc);

  this->Configure(init_state);
}
//___________________________________________________________________________
void GEVGDriver::Configure(const InitialState & init_state)
{
  string mesgh = "Configuring a GEVGDriver object for init-state = ";
  LOG("GEVGDriver", pNOTICE) 
    << utils::print::PrintFramedMesg(mesgh + init_state.AsString(), 0, '*');

  this -> BuildInitialState            (init_state);
  this -> BuildGeneratorList           ();
  this -> BuildInteractionGeneratorMap ();
  this -> BuildInteractionSelector     ();

  LOG("GEVGDriver", pNOTICE) << "Done configuring GEVGDriver! \n";
}
//___________________________________________________________________________
void GEVGDriver::BuildInitialState(const InitialState & init_state)
{
  LOG("GEVGDriver", pNOTICE) << "Setting the initial state...";

  if(fInitState) delete fInitState;
  fInitState = new InitialState(init_state);

  this->AssertIsValidInitState();
}
//___________________________________________________________________________
void GEVGDriver::BuildGeneratorList(void)
{
//! figure out which list of event generators to use from the $GEVGL
//! environmental variable (use "Default") if the variable is not set.

  LOG("GEVGDriver", pNOTICE) << "Building the event generator list...";

  string evgl = (gSystem->Getenv("GEVGL") ?
                               gSystem->Getenv("GEVGL") : "Default");
  LOG("GEVGDriver", pNOTICE)
                       << "Specified Event Generator List = " << evgl;

  EventGeneratorListAssembler evglist_assembler(evgl.c_str());
  fEvGenList = evglist_assembler.AssembleGeneratorList();
}
//___________________________________________________________________________
void GEVGDriver::BuildInteractionGeneratorMap(void)
{
//! figure out which list of event generators to use from the $GEVGL
//! environmental variable (use "Default") if the variable is not set.

  LOG("GEVGDriver", pNOTICE)
         << "Building the interaction -> generator associations...";

  fIntGenMap = new InteractionGeneratorMap;
  fIntGenMap->UseGeneratorList(fEvGenList);
  fIntGenMap->BuildMap(*fInitState);

  string mesgh = "Interaction -> Generator assignments for Initial State: ";

  LOG("GEVGDriver", pNOTICE)
    << utils::print::PrintFramedMesg(mesgh + fInitState->AsString(), 0, '-')
    << *fIntGenMap;
}
//___________________________________________________________________________
void GEVGDriver::BuildInteractionSelector(void)
{
  LOG("GEVGDriver", pNOTICE) << "Building the interaction selector...";

  if(fIntSelector) delete fIntSelector;
  fIntSelector = new PhysInteractionSelector("Default");
}
//___________________________________________________________________________
EventRecord * GEVGDriver::GenerateEvent(const TLorentzVector & nu4p)
{
  //-- Build initial state information from inputs
  LOG("GEVGDriver", pINFO) << "Creating the initial state";
  InitialState init_state(*fInitState);
  init_state.SetProbeP4(nu4p);

  //-- Select the interaction to be generated (amongst the entries of the
  //   InteractionList assembled by the EventGenerators) and bootstrap the
  //   event record
  LOG("GEVGDriver", pINFO)
             << "Selecting an Interaction & Bootstraping the EventRecord";
  fCurrentRecord = fIntSelector->SelectInteraction(fIntGenMap, nu4p);

  assert(fCurrentRecord); // abort if no interaction could be selected!

  //-- Get a ptr to the interaction summary
  LOG("GEVGDriver", pDEBUG) << "Getting the selected interaction";
  Interaction * interaction = fCurrentRecord->Summary();

  //-- Find the appropriate concrete EventGeneratorI implementation
  //   for generating this event.
  //
  //   The right EventGeneratorI will be selecting by iterating over the
  //   entries of the EventGeneratorList and compare the interaction
  //   against the ValidityContext declared by each EventGeneratorI
  //
  //   (note: use of the 'Chain of Responsibility' Design Pattern)

  LOG("GEVGDriver", pINFO) << "Finding an appropriate EventGenerator";

  const EventGeneratorI * evgen = fIntGenMap->FindGenerator(interaction);
  assert(evgen);

  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  rtinfo->UpdateRunningThread(evgen);

  //-- Generate the selected event
  //
  //   The selected EventGeneratorI subclass will start processing the
  //   event record (by sequentially asking each entry in its list of
  //   EventRecordVisitorI subclasses to visit and process the record).
  //   Most of the actual event generation takes place in this step.
  //
  //   (note: use of the 'Visitor' Design Pattern)

  string mesg = "Requesting from event generation thread: " + 
         evgen->Id().Key() + " to generate the selected interaction";

  LOG("GEVGDriver", pNOTICE) 
         << utils::print::PrintFramedMesg(mesg,1,'=');

  evgen->ProcessEventRecord(fCurrentRecord);

  //-- Check whether the generated event flags. The default behaviour is
  //   to filter out unphysical events and enter in recursive mode to 
  //   regenerate them.
  //   If an unphysical event mask has been set, apply it to the event
  //   flags and allow the requested classes of unphysical events to be 
  //   returned

  bool unphys = fCurrentRecord->IsUnphysical();
  if(unphys) {
     LOG("GEVGDriver", pWARN) << "I generated an unphysical event!";
  }

  TBits evflags = *(fCurrentRecord->EventFlags());
  TBits mask    = *(fFiltUnphysMask);

  TBits flags   = evflags & (~mask);
  bool  filter  = (flags.CountBits()>0);

  if(!filter) {
     LOG("GEVGDriver", pNOTICE) << "Returning the current event!";
     fNRecLevel = 0;     
     return fCurrentRecord; // The client 'adopts' the event record

  } else {
     LOG("GEVGDriver", pWARN) << "I am filtering out the current event!";
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
        exit(1);
     }
  }
}
//___________________________________________________________________________
const InteractionList * GEVGDriver::Interactions(void) const
{
// Returns the list of all interactions that can be generated by this driver

  if(!fIntGenMap) {
    LOG("GEVGDriver", pWARN)
       << "Interaction->Generator Map has not being built yet!";
    return 0;
  }

  const InteractionList & ilst = fIntGenMap->GetInteractionList();
  return &ilst;
}
//___________________________________________________________________________
const EventGeneratorI * GEVGDriver::FindGenerator(
                                  const Interaction * interaction) const
{
  if(!interaction) {
    LOG("GEVGDriver", pWARN) << "Null interaction!!";
    return 0;
  }
  if(!fIntGenMap) {
    LOG("GEVGDriver", pWARN)
       << "Interaction->Generator Map has not being built yet!";
    return 0;
  }

  const EventGeneratorI * evgen = fIntGenMap->FindGenerator(interaction);
  return evgen;
}
//___________________________________________________________________________
double GEVGDriver::XSecSum(const TLorentzVector & nup4)
{
// Computes the sum of the cross sections for all the interactions that can
// be simulated for the given initial state and for the input neutrino energy
//
  LOG("GEVGDriver", pDEBUG) << "Computing the cross section sum";

  double xsec_sum = 0;

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = XSecSplineList::Instance();

  // Get the list of all interactions that can be generated by this driver
  const InteractionList & ilst = fIntGenMap->GetInteractionList();

  // Loop over all interactions & compute cross sections
  InteractionList::const_iterator intliter;
  for(intliter = ilst.begin(); intliter != ilst.end(); ++intliter) {

     // get current interaction
     Interaction * interaction = new Interaction(**intliter);
     interaction->InitStatePtr()->SetProbeP4(nup4);

     string code = interaction->AsString();
     SLOG("GEVGDriver", pDEBUG)
             << "Compute cross section for interaction: \n" << code;

     // get corresponding cross section algorithm
     const XSecAlgorithmI * xsec_alg =
               fIntGenMap->FindGenerator(interaction)->CrossSectionAlg();
     assert(xsec_alg);

     // compute (or evaluate) the cross section
     double xsec = 0;
     bool spline_exists = xssl->SplineExists(xsec_alg, interaction);
     if (spline_exists && fUseSplines) {
        double E = nup4.Energy();
        xsec = xssl->GetSpline(xsec_alg,interaction)->Evaluate(E);
     } else
        xsec = xsec_alg->Integral(interaction);

     // sum-up and report
     xsec_sum += xsec;
     LOG("GEVGDriver", pDEBUG)
            << "\nInteraction   = " << code
            << "\nCross Section "
            << (fUseSplines ? "*interpolated*" : "*computed*")
            << " = " << (xsec/units::cm2) << " cm2";

     delete interaction;
  } // loop over event generators

  PDGLibrary * pdglib = PDGLibrary::Instance();
  LOG("GEVGDriver", pINFO)
    << "SumXSec("
    << pdglib->Find(fInitState->ProbePdg())->GetName() << "+"
    << pdglib->Find(fInitState->Tgt().Pdg())->GetName() << "->X, "
    << "E = " << nup4.Energy() << " GeV)"
    << (fUseSplines ? "*interpolated*" : "*computed*")
    << " = " << (xsec_sum/units::cm2) << " cm2";

  return xsec_sum;
}
//___________________________________________________________________________
void GEVGDriver::CreateXSecSumSpline(
                               int nk, double Emin, double Emax, bool inlogE)
{
// This method creates a spline with the *total* cross section vs E (or logE)
// for the initial state that this driver was configured with.
// This spline is used, for example, by the GMCJDriver to select a target
// material out of all the materials in a detector geometry (summing the
// cross sections again and again proved to be expensive...)

  LOG("GEVGDriver", pNOTICE)
     << "Creating spline (sum-xsec = f(" << ((inlogE) ? "logE" : "E")
     << ") in E = [" << Emin << ", " << Emax << "] using " << nk << " knots";

  if(!fUseSplines) {
     LOG("GEVGDriver", pFATAL) << "You haven't loaded any splines!! ";
  }
  assert(fUseSplines);
  assert(Emin<Emax && Emin>0 && nk>2);

  double logEmin=0, logEmax=0, dE=0;

  double * E    = new double[nk];
  double * xsec = new double[nk];

  if(inlogE) {
    logEmin = TMath::Log(Emin);
    logEmax = TMath::Log(Emax);
    dE = (logEmax-logEmin)/(nk-1);
  } else {
    dE = (Emax-Emin)/(nk-1);
  }

  TLorentzVector p4(0,0,0,0);

  for(int i=0; i<nk; i++) {
    double e = (inlogE) ? TMath::Exp(logEmin + i*dE) : Emin + i*dE;
    p4.SetPxPyPzE(0.,0.,e,e);
    double xs = this->XSecSum(p4);

    E[i]    = e;
    xsec[i] = xs;
  }
  if (fXSecSumSpl) delete fXSecSumSpl;
  fXSecSumSpl = new Spline(nk, E, xsec);
}
//___________________________________________________________________________
const Spline * GEVGDriver::XSecSpline(const Interaction * interaction) const
{
// Returns the cross section spline for the input interaction as was
// computed from the cross section model associated with that interaction.

  if (!fUseSplines) return 0;

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xssl = XSecSplineList::Instance();

  // get corresponding cross section algorithm for the input interaction
  const XSecAlgorithmI * xsec_alg =
               fIntGenMap->FindGenerator(interaction)->CrossSectionAlg();
  assert(xsec_alg);

  const Spline * spl = xssl->GetSpline(xsec_alg,interaction);
  return spl;
}
//___________________________________________________________________________
void GEVGDriver::UseSplines(void)
{
// Instructs the driver to use cross section splines rather than computing
// cross sections by integrating the differential cross section model which
// can be very time-consuming.
// **Note**
// -- If you called GEVGDriver::CreateSplines() already the driver would
//    a) assume that you want to use them and b) would be assured that it 
//    has all the splines it needs, so you do not need to call this method.
// -- If you populated the XSecSplineList in another way without this driver
//    knowing about it, eg from an external XML file, do call this method 
//    to let the driver know that you want to use the splines. However, note
//    that the driver would **explicitly check** that you have loaded all the
//    splines it needs. If not, then its fiery personality will take over and
//    it will refuse your request, reverting back to not using splines.

  fUseSplines = true;

  // Get the list of spline objects
  // Should have been constructed at the job initialization
  XSecSplineList * xsl = XSecSplineList::Instance();

  // If the user wants to use splines, make sure that all the splines needed
  // have been computed or loaded
  if(fUseSplines) {

    // Get the list of all interactions that can be generated by this driver
    const InteractionList & ilst = fIntGenMap->GetInteractionList();

    // Loop over all interactions & check that all splines have been loaded
    InteractionList::const_iterator intliter;
    for(intliter = ilst.begin(); intliter != ilst.end(); ++intliter) {

       // current interaction
       Interaction * interaction = *intliter;

       // corresponding cross section algorithm
       const XSecAlgorithmI * xsec_alg =
                 fIntGenMap->FindGenerator(interaction)->CrossSectionAlg();
       assert(xsec_alg);

       // spline exists in spline list?
       bool spl_exists = xsl->SplineExists(xsec_alg, interaction);

       // update the 'use splines' flag
       fUseSplines = fUseSplines && spl_exists;

       if(!spl_exists) {
          LOG("GEVGDriver", pWARN)
             << "*** At least a spline doesn't exist. "
                               << "Reverting back to not using splines";
          return;
       }
     } // loop over interaction list
  }//use-splines?
}
//___________________________________________________________________________
void GEVGDriver::CreateSplines(int nknots, double emax, bool useLogE)
{
// Creates all the cross section splines that are needed by this driver.
// It will check for pre-loaded splines and it will skip the creation of the
// splines it already finds loaded.

  LOG("GEVGDriver", pINFO)
       << "Creating (missing) splines with [UseLogE: "
                                             << ((useLogE) ? "ON]" : "OFF]");
  // Get the list of spline objects
  XSecSplineList * xsl = XSecSplineList::Instance();
  xsl->SetLogE(useLogE);

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  // loop over all EventGenerator objects used in the current job
  for(evgliter = fEvGenList->begin();
                               evgliter != fEvGenList->end(); ++evgliter) {
     // current event generator
     const EventGeneratorI * evgen = *evgliter;
     LOG("GEVGDriver", pNOTICE)
            << "Querying [" << evgen->Id().Key()
                                            << "] for its InteractionList";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state
     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(*fInitState);
     if(!ilst) continue;

     // total cross section algorithm used by the current EventGenerator
     const XSecAlgorithmI * alg = evgen->CrossSectionAlg();

     // get the energy range of the spline from the EventGenerator
     // validity context
     double Emin = TMath::Max(0.01,evgen->ValidityContext().Emin());
     double Emax = evgen->ValidityContext().Emax();

     // if the user set a maximum energy, create the spline up to this
     // energy - otherwise use the upper limit of the validity range of
     // the current generator
     if(emax>0) {
       if(emax>Emax) {
	 LOG("GEVGDriver", pWARN) 
           << "Refusing to exceed validity range: Emax = " << Emax;
       }
       emax = TMath::Min(emax,Emax); // don't exceed validity range
     } else emax = Emax;

     assert(emax>Emin);

     // number of knots: use specified number. If not set, use 15 knots
     // per decade. Don't use less than 30 knots.
     if(nknots<0) {
       nknots = (int) (15 * TMath::Log10(emax-Emin));
     }
     nknots = TMath::Max(nknots,30);

     // loop over all interactions that can be generated and ask the
     // appropriate cross section algorithm to compute its cross section
     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter) {

         // current interaction
         Interaction * interaction = *intliter;
         string code = interaction->AsString();

         SLOG("GEVGDriver", pNOTICE) << "Need xsec spline for " << code;

         // only create the spline if it does not already exists
         bool spl_exists = xsl->SplineExists(alg, interaction);
         if(!spl_exists) {
             SLOG("GEVGDriver", pNOTICE) 
               << "The spline wasn't loaded at initialization. "
               << "I can build it now but it might take a while..."; 
             SLOG("GEVGDriver", pNOTICE) 
               << "(Consider recycling the splines! "
               << "See the GENIE FAQ/HOWTO at http://www.genie-mc.org)";
             xsl->CreateSpline(alg, interaction, nknots, Emin, emax);
         } else {
             SLOG("GEVGDriver", pNOTICE) << "Spline was found";
         }
     } // loop over interaction that can be generated by this generator
     delete ilst;
     ilst = 0;
  } // loop over event generators

  LOG("GEVGDriver", pINFO) << *xsl; // print list of splines

  fUseSplines = true;
}
//___________________________________________________________________________
Range1D_t GEVGDriver::ValidEnergyRange(void) const
{
// loops over all loaded event generation threads, queries for the energy
// range at their validity context and builds the valid energy range for
// this driver

  Range1D_t E;
  E.min =  9999;
  E.max = -9999;

  EventGeneratorList::const_iterator evgliter; // event generator list iter

  // loop over all EventGenerator objects used in the current job
  for(evgliter = fEvGenList->begin();
                               evgliter != fEvGenList->end(); ++evgliter) {
     // current event generator
     const EventGeneratorI * evgen = *evgliter;

     // Emin, Emax as declared in current generator's validity context
     double Emin = TMath::Max(0.01,evgen->ValidityContext().Emin());
     double Emax = evgen->ValidityContext().Emax();

     // combined Emin, Emax
     E.min = TMath::Min(E.min, Emin);
     E.max = TMath::Max(E.max, Emax);
  }
  assert(E.min<E.max && E.min>=0);
  return E;
}
//___________________________________________________________________________
void GEVGDriver::AssertIsValidInitState(void) const
{
  assert(fInitState);

  int nu_pdgc = fInitState->ProbePdg();

  bool isnu = pdg::IsNeutrino(nu_pdgc) || pdg::IsAntiNeutrino(nu_pdgc);
  assert(isnu);
}
//___________________________________________________________________________
void GEVGDriver::Print(ostream & stream) const
{
  stream
    << "\n\n *********************** GEVGDriver ***************************";

  int nupdg  = fInitState->ProbePdg();
  int tgtpdg = fInitState->Tgt().Pdg();

  stream << "\n  |---o Neutrino PDG-code .........: " << nupdg;
  stream << "\n  |---o Nuclear Target PDG-code ...: " << tgtpdg;

  stream << "\n  |---o Using cross section splines is turned "
                                << utils::print::BoolAsIOString(fUseSplines);
  stream << "\n  |---o Unphysical event filter mask ("
         << GHepFlags::NFlags() << "->0) = " << *fFiltUnphysMask;

  stream << "\n *********************************************************\n";
}
//___________________________________________________________________________

