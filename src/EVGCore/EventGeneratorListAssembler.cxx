//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - January 25, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventGeneratorListAssembler.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "Messenger/Messenger.h"
#include "Utils/PrintUtils.h"

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
  SLOG("EvGenListAssembler", pNOTICE) 
            << utils::print::PrintFramedMesg(
                           "Loading requested Event Generators", 0, '-');

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

  SLOG("EvGenListAssembler", pNOTICE) 
        << "Loading generator: " << alg << "/" << conf;

  AlgFactory * algf = AlgFactory::Instance();

  const EventGeneratorI * evgen =
      dynamic_cast<const EventGeneratorI *> (algf->GetAlgorithm(alg,conf));
  assert(evgen);

  SLOG("EvGenListAssembler", pNOTICE) << "\n";

  return evgen;
}
//___________________________________________________________________________

