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

#include <sstream>

#include <TMCParticle6.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EventGenerator.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "GHEP/GHepVirtualListFolder.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepRecordHistory.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
EventGenerator::EventGenerator() :
EventGeneratorI()
{
  fName = "genie::EventGenerator";
}
//___________________________________________________________________________
EventGenerator::EventGenerator(const char * param_set) :
EventGeneratorI(param_set)
{
  fName = "genie::EventGenerator";

  this->FindConfig();

  this->InstantiateValidityContext();
}
//___________________________________________________________________________
EventGenerator::~EventGenerator()
{

}
//___________________________________________________________________________
void EventGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  SLOG("EventGenerator", pINFO) << "Generating Event:";

  fConfig->AssertExistence("n-generator-steps");
  int nsteps = fConfig->GetInt("n-generator-steps");

  if(nsteps == 0) {
    LOG("EventGenerator", pWARN) 
         << "EventGenerator configuration declares null visitor list!";
  }

  //-- Clear previous virtual list folder
  LOG("EventGenerator", pINFO) << "Clearing the GHepVirtualListFolder";
  GHepVirtualListFolder * vlfolder = GHepVirtualListFolder::Instance();
  vlfolder->Clear();

  //-- Create a history buffer in case I need to step back
  GHepRecordHistory rh;

  //-- Loop over the event record processing steps
  for(int istep = 0; istep < nsteps; istep++) {

     const EventRecordVisitorI * visitor = this->ProcessingStep(istep);

     bool ffwd = event_rec->FastForwardEnabled();
     if(!ffwd) {
         visitor->ProcessEventRecord(event_rec);
         rh.AddSnapshot(istep, event_rec);
     } else {
       LOG("EventGenerator", pINFO)
           << "Fast Forward flag was set - Skipping processing step!";
     }
  }

  LOG("EventGenerator", pINFO)
         << "The EventRecord was visited by all EventRecordVisitors\n";
}
//___________________________________________________________________________
const EventRecordVisitorI * EventGenerator::ProcessingStep(int istep) const
{
  ostringstream alg_key, config_key;

  alg_key    << "generator-step-" << istep << "-alg";
  config_key << "generator-step-" << istep << "-conf";

  string alg    = fConfig->GetString( alg_key.str()    );
  string config = fConfig->GetString( config_key.str() );

  SLOG("EventGenerator", pINFO)
      << "\n\n---------- Getting EventRecordVisitorI algorithm: "
                                  << alg << "/" << config << " ----------\n";

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg, config);

  const EventRecordVisitorI * visitor =
                        dynamic_cast<const EventRecordVisitorI *> (algbase);

  return visitor;
}
//___________________________________________________________________________
const InteractionListGeneratorI * EventGenerator::IntListGenerator(void) const
{
  const Algorithm * algbase = this->SubAlg(
                            "interaction-list-alg", "interaction-list-conf");
  const InteractionListGeneratorI * intlistgen =
                   dynamic_cast<const InteractionListGeneratorI *> (algbase);
  assert(intlistgen);

  return intlistgen;
}
//___________________________________________________________________________
const XSecAlgorithmI * EventGenerator::CrossSectionAlg(void) const
{
  const Algorithm * algbase = this->SubAlg(
                                  "cross-section-alg", "cross-section-conf");
  const XSecAlgorithmI * xsecalg =
                              dynamic_cast<const XSecAlgorithmI *> (algbase);
  assert(xsecalg);

  return xsecalg;
}
//___________________________________________________________________________
void EventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);

  this->InstantiateValidityContext();
}
//___________________________________________________________________________
void EventGenerator::Configure(string param_set)
{
  Algorithm::Configure(param_set);

  this->InstantiateValidityContext();
}
//___________________________________________________________________________
const GVldContext & EventGenerator::ValidityContext(void) const
{
  return *fVldContext;
}
//___________________________________________________________________________
void EventGenerator::InstantiateValidityContext(void)
{
  fVldContext = new GVldContext;

  assert( fConfig->Exists("vld-context") );

  string encoded_vld_context = fConfig->GetString("vld-context");

  fVldContext->Decode( encoded_vld_context );
}
//___________________________________________________________________________
