//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/EventGeneratorListAssembler.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/RunOpt.h"

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

  if ( RunOpt::Instance() -> EventGeneratorList() == "Default" ) {
    AddTopRegistry( AlgConfigPool::Instance() -> TuneGeneratorList(), false ) ;
    
    SLOG("EvGenListAssembler", pNOTICE) 
      << "** Using Tune Generator List: " ;

  }
  

  EventGeneratorList * evgl = new EventGeneratorList;
  
//  if (!fConfig) {
//    SLOG("EvGenListAssembler", pFATAL)
//      << "Cannot instantiate EventGeneratorList with no config.";
//    gAbortingInErr = true;
//    exit(-1);
//  }

  int nproc = GetConfig().GetInt("NGenerators");
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

