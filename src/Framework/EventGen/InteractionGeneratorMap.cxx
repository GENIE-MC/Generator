//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <iomanip>

#include <TMath.h>

#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/InteractionGeneratorMap.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

using std::setw;
using std::setfill;
using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const InteractionGeneratorMap & intl)
 {
   intl.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
InteractionGeneratorMap::InteractionGeneratorMap() :
map<string, const EventGeneratorI *> ()
{
  this->Init();
}
//___________________________________________________________________________
InteractionGeneratorMap::InteractionGeneratorMap(
                                     const InteractionGeneratorMap & igmap) :
map<string, const EventGeneratorI *> ()
{
  this->Copy(igmap);
}
//___________________________________________________________________________
InteractionGeneratorMap::~InteractionGeneratorMap()
{
  this->CleanUp();
}
//___________________________________________________________________________
void InteractionGeneratorMap::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void InteractionGeneratorMap::Init(void)
{
  fEventGeneratorList = 0;

  fInitState       = new InitialState;
  fInteractionList = new InteractionList;
}
//___________________________________________________________________________
void InteractionGeneratorMap::CleanUp(void)
{
  delete fInitState;
  delete fInteractionList;

  this->clear();
}
//___________________________________________________________________________
void InteractionGeneratorMap::Copy(const InteractionGeneratorMap & xsmap)
{
  fEventGeneratorList = xsmap.fEventGeneratorList;

  fInitState       -> Copy (*xsmap.fInitState);
  fInteractionList -> Copy (*xsmap.fInteractionList);

  this->clear();

  InteractionGeneratorMap::const_iterator iter;

  for(iter = xsmap.begin(); iter != xsmap.end(); ++iter) {
    string code = iter->first;
    const EventGeneratorI * evg = iter->second;

    this->insert(map<string, const EventGeneratorI *>::value_type(code,evg));
  }
}
//___________________________________________________________________________
void InteractionGeneratorMap::UseGeneratorList(const EventGeneratorList * l)
{
  fEventGeneratorList = l;
}
//___________________________________________________________________________
void InteractionGeneratorMap::BuildMap(const InitialState & init_state)
{
  SLOG("IntGenMap", pDEBUG)
               << "Building 'interaction' -> 'generator' associations";
  SLOG("IntGenMap", pNOTICE)
         << "Using all simulated interactions for init-state: "
                                              << init_state.AsString();
  if(!fEventGeneratorList) {
     LOG("IntGenMap", pWARN) << "No EventGeneratorList was loaded!!";
     return;
  }

  fInitState->Copy(init_state);

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  // loop over all EventGenerator objects used in the current job
  for(evgliter = fEventGeneratorList->begin();
                       evgliter != fEventGeneratorList->end(); ++evgliter) {
     // current EventGenerator
     const EventGeneratorI * evgen = *evgliter;
     assert(evgen);

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state
     SLOG("IntGenMap", pNOTICE)
        << "Querying [" << evgen->Id().Key() << "] for its InteractionList";

     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(init_state);

     // no point to go on if the list is NULL - continue to next iteration
     if(!ilst) continue;

     // append the new InteractionList to the local copy
     fInteractionList->Append(*ilst);

     // loop over all interaction that can be genererated by the current
     // EventGenerator and link all of them to iy
     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter)
     {
        // current interaction
        Interaction * interaction = *intliter;
        string code = interaction->AsString();

        SLOG("IntGenMap", pDEBUG)
              << "\nLinking: " << code << " --> to: " << evgen->Id().Key();
        this->insert(
             map<string, const EventGeneratorI *>::value_type(code,evgen));
     } // loop over interactions
     delete ilst;
     ilst = 0;
  } // loop over event generators
}
//___________________________________________________________________________
const EventGeneratorI * InteractionGeneratorMap::FindGenerator(
                                      const Interaction * interaction) const
{
  if(!interaction) {
    LOG("IntGenMap", pWARN) << "Null interaction!!";
    return 0;
  }
  string code = interaction->AsString();
  InteractionGeneratorMap::const_iterator evgiter = this->find(code);
  if(evgiter == this->end()) {
    LOG("IntGenMap", pWARN)
             << "No EventGeneratorI was found for interaction: \n" << code;
    return 0;
  }
  const EventGeneratorI * evg = evgiter->second;
  return evg;
}
//___________________________________________________________________________
const InteractionList & InteractionGeneratorMap::GetInteractionList(void) const
{
  return *fInteractionList;
}
//___________________________________________________________________________
void InteractionGeneratorMap::Print(ostream & stream) const
{
  stream << endl;

  InteractionGeneratorMap::const_iterator iter;

  unsigned int maxlen = 0;
  for(iter = this->begin(); iter != this->end(); ++iter) {
    string icode = iter->first;
    unsigned int isz = (unsigned int) icode.size();
    maxlen=TMath::Max(maxlen,isz);
  }

  for(iter = this->begin(); iter != this->end(); ++iter) {
    const EventGeneratorI * evg = iter->second;
    string intstr = iter->first;
    string evgstr = (evg) ? evg->Id().Key() : "** NULL EVENT GENERATOR **";

    stream << setfill(' ') << setw(maxlen) 
           << intstr << " --> " << evgstr << endl;
  }
}
//___________________________________________________________________________
InteractionGeneratorMap & InteractionGeneratorMap::operator = (
                                       const InteractionGeneratorMap & igmap)
{
  this->Copy(igmap);
  return (*this);
}
//___________________________________________________________________________


