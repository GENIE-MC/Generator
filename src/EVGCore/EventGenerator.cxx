//____________________________________________________________________________
/*!

\class   genie::EventGenerator

\brief   Encapsulates a full ordered list of (is the aggregate of) concrete
         EventGeneratorI implementations that must act on the EventRecord
         to generate an event. Each of these implementations corresponds to
         a single processing step.

         Is a concrete implementation of the EventGeneratorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>

#include <TStopwatch.h>
#include <TMCParticle6.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventGenerator.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/GVldContext.h"
#include "GHEP/GHepVirtualListFolder.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepRecordHistory.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace genie::exceptions;

//___________________________________________________________________________
EventGenerator::EventGenerator() :
EventGeneratorI("genie::EventGenerator")
{

}
//___________________________________________________________________________
EventGenerator::EventGenerator(string config) :
EventGeneratorI("genie::EventGenerator", config)
{

}
//___________________________________________________________________________
EventGenerator::~EventGenerator()
{
  delete fWatch;

  if(fEVGModuleVec) delete fEVGModuleVec;
  if(fEVGTime)      delete fEVGTime;
  if(fVldContext)   delete fVldContext;
}
//___________________________________________________________________________
void EventGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  SLOG("EventGenerator", pNOTICE) << "Generating Event:";

  //-- Clear previous virtual list folder
  LOG("EventGenerator", pNOTICE) << "Clearing the GHepVirtualListFolder";
  GHepVirtualListFolder * vlfolder = GHepVirtualListFolder::Instance();
  vlfolder->Clear();

  //-- Create a history buffer in case I need to step back
  GHepRecordHistory rh;
  rh.AddSnapshot(-1, event_rec);

  //-- Initialize evg thread control flags
  bool ffwd = false;
  unsigned int nexceptions = 0;

  //-- Reset stop-watch
  fWatch->Reset();

  //-- Loop over the event record processing steps
  int istep=0;
  vector<const EventRecordVisitorI *>::const_iterator miter;

  for(miter = fEVGModuleVec->begin();
                               miter != fEVGModuleVec->end(); ++miter){

    const EventRecordVisitorI * visitor = *miter; // generation module

    SLOG("EventGenerator", pNOTICE)
      << "\n\n---------- Running EventRecordVisitorI: "
                             << visitor->Id().Key() << " ----------\n";
    if(ffwd) {
      LOG("EventGenerator", pNOTICE)
           << "Fast Forward flag was set - Skipping processing step!";
      continue;
    }
    try
    {
      fWatch->Start();
      visitor->ProcessEventRecord(event_rec);
      fWatch->Stop();

      rh.AddSnapshot(istep, event_rec);

      (*fEVGTime)[istep] = fWatch->CpuTime(); // sec
    }
    catch (EVGThreadException exception)
    {
      LOG("EventGenerator", pNOTICE)
           << "An exception was thrown and caught by EventGenerator!";
      LOG("EventGenerator", pNOTICE) << exception;

      nexceptions++;
      if ( nexceptions > kMaxEVGThreadExceptions ) {
         LOG("EventGenerator", pFATAL)
           << "Caught max allowed number (" << kMaxEVGThreadExceptions
                           << ") of EVGThreadExceptions/thread. Aborting";
         abort();
      }

      // make sure we are not asked to go at both directions...
      assert( !(exception.FastForward() && exception.StepBack()) );

      ffwd = exception.FastForward();

      if(exception.StepBack()) {

         // get return step (if return_step > current_step just ignore it)
         if(exception.ReturnStep() >= 0 && exception.ReturnStep() <= istep) {
           istep = exception.ReturnStep();
           istep--;

           // restore the event record as it was just before the processing
           // step we are about to return to
           event_rec->ResetRecord();
           GHepRecord * snapshot = rh[istep];
           rh.PurgeRecentHistory(istep+1);
           event_rec->Copy(*snapshot);
         } // valid-return-step
      } // step-back
    } // catch exception

    istep++;
  }

  LOG("EventGenerator", pNOTICE)
         << "The EventRecord was visited by all EventRecordVisitors\n";

  LOG("EventGenerator", pINFO) << "** Event generation timing info **";
  istep=0;
  for(miter = fEVGModuleVec->begin();
                               miter != fEVGModuleVec->end(); ++miter){
    const EventRecordVisitorI * visitor = *miter;

    SLOG("EventGenerator", pINFO)
      << "EventRecordVisitorI: " << visitor->Id().Key() << " -> "
                                          << (*fEVGTime)[istep++] << " sec";
  }
}
//___________________________________________________________________________
const InteractionListGeneratorI * EventGenerator::IntListGenerator(void) const
{
  return fIntListGen;
}
//___________________________________________________________________________
const XSecAlgorithmI * EventGenerator::CrossSectionAlg(void) const
{
  return fXSecModel;
}
//___________________________________________________________________________
void EventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);

  this->Init();
  this->LoadVldContext();
  this->LoadEVGModules();
  this->LoadIntSelAlg();
}
//___________________________________________________________________________
void EventGenerator::Configure(string param_set)
{
  Algorithm::Configure(param_set);

  this->Init();
  this->LoadVldContext();
  this->LoadEVGModules();
  this->LoadIntSelAlg();
}
//___________________________________________________________________________
const GVldContext & EventGenerator::ValidityContext(void) const
{
  return *fVldContext;
}
//___________________________________________________________________________
void EventGenerator::Init(void)
{
  fWatch        = new TStopwatch;
  fVldContext   = 0;
  fEVGModuleVec = 0;
  fEVGTime      = 0;
  fXSecModel    = 0;
  fIntListGen   = 0;
}
//___________________________________________________________________________
void EventGenerator::LoadVldContext(void)
{
  fVldContext = new GVldContext;

  assert( fConfig->Exists("vld-context") );
  string encoded_vld_context = fConfig->GetString("vld-context");

  fVldContext->Decode( encoded_vld_context );
}
//___________________________________________________________________________
void EventGenerator::LoadEVGModules(void)
{
  fConfig->AssertExistence("n-generator-steps");
  int nsteps = fConfig->GetInt("n-generator-steps");

  if(nsteps == 0) {
    LOG("EventGenerator", pFATAL)
         << "EventGenerator configuration declares null visitor list!";
  }
  assert(nsteps>0);

  fEVGModuleVec = new vector<const EventRecordVisitorI *> (nsteps);
  fEVGTime      = new vector<double>(nsteps);

  for(int istep = 0; istep < nsteps; istep++) {

     ostringstream alg_key, config_key;

     alg_key    << "generator-step-" << istep << "-alg";
     config_key << "generator-step-" << istep << "-conf";

     string alg    = fConfig->GetString( alg_key.str()    );
     string config = fConfig->GetString( config_key.str() );
     SLOG("EventGenerator", pNOTICE)
           << "Processing step: " << istep
               << "-> loading EventRecordVisitorI: " << alg << "/" << config;

     const EventRecordVisitorI * visitor =
                dynamic_cast<const EventRecordVisitorI *> (
                              this->SubAlg(alg_key.str(), config_key.str()));
     (*fEVGModuleVec)[istep] = visitor;
     (*fEVGTime)[istep]      = 0;
  }
}
//___________________________________________________________________________
void EventGenerator::LoadIntSelAlg(void)
{
  fIntListGen = dynamic_cast<const InteractionListGeneratorI *> (
              this->SubAlg("interaction-list-alg", "interaction-list-conf"));
  assert(fIntListGen);

  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                    this->SubAlg("cross-section-alg", "cross-section-conf"));
  assert(fXSecModel);
}
//___________________________________________________________________________


