//____________________________________________________________________________
/*!

\class   genie::EventGeneratorListAssembler

\brief   Assembles a list of all the EventGeneratorI subclasses that can be
         employed during a neutrino event generation job.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created January 25, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventGeneratorListAssembler.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
EventGeneratorListAssembler::EventGeneratorListAssembler() :
Algorithm("genie::EventGeneratorListAssembler")
{

}
//___________________________________________________________________________
EventGeneratorListAssembler::EventGeneratorListAssembler(string config) :
Algorithm("genie::EventGeneratorListAssembler", config)
{

}
//___________________________________________________________________________
EventGeneratorListAssembler::~EventGeneratorListAssembler()
{

}
//___________________________________________________________________________
EventGeneratorList * EventGeneratorListAssembler::AssembleGeneratorList()
{
  SLOG("EvGenListAssembler", pNOTICE) << "----------------------------------";
  SLOG("EvGenListAssembler", pNOTICE) << "Loading requested Event Generators";
  SLOG("EvGenListAssembler", pNOTICE) << "----------------------------------";

  EventGeneratorList * evgl = new EventGeneratorList;

  int nproc = fConfig->GetIntDef("n-processes", 0);
  assert(nproc > 0);

  //-- Loop over the event generators for all requested processes
  for(int ip = 0; ip < nproc; ip++) {
    const EventGeneratorI * evgen = this->LoadGenerator(ip);
    evgl->push_back(evgen);
  }

  return evgl;
}
//___________________________________________________________________________
const EventGeneratorI * EventGeneratorListAssembler::LoadGenerator(int ip)
{
  ostringstream alg_key, config_key;

  alg_key    << "event-generator-id-" << ip << "-alg-name";
  config_key << "event-generator-id-" << ip << "-param-set";

  string alg  = fConfig->GetString( alg_key.str()    );
  string conf = fConfig->GetString( config_key.str() );

  SLOG("EvGenListAssembler", pNOTICE) << "Loading: " << alg << "/" << conf;

  AlgFactory * algf = AlgFactory::Instance();

  const EventGeneratorI * evgen =
      dynamic_cast<const EventGeneratorI *> (algf->GetAlgorithm(alg,conf));
  assert(evgen);

  SLOG("EvGenListAssembler", pNOTICE) << "\n";

  return evgen;
}
//___________________________________________________________________________

