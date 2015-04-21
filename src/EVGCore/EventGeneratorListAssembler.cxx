//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - January 25, 2004

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

  int nproc = fConfig->GetIntDef("NGenerators", 0);
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
  ostringstream alg_key;
  alg_key << "Generator-" << ip;

  const EventGeneratorI * evgen =
      dynamic_cast<const EventGeneratorI *> (this->SubAlg(alg_key.str()));
  assert(evgen);

  SLOG("EvGenListAssembler", pNOTICE) 
        << "** Loaded generator: " << evgen->Id().Key();

  return evgen;
}
//___________________________________________________________________________

