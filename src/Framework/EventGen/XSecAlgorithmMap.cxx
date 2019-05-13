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

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/XSecAlgorithmMap.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const XSecAlgorithmMap & intl)
 {
   intl.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
XSecAlgorithmMap::XSecAlgorithmMap() :
map<string, const XSecAlgorithmI *> ()
{
  this->Init();
}
//___________________________________________________________________________
XSecAlgorithmMap::XSecAlgorithmMap(const XSecAlgorithmMap & xsmap) :
map<string, const XSecAlgorithmI *> ()
{
  this->Copy(xsmap);
}
//___________________________________________________________________________
XSecAlgorithmMap::~XSecAlgorithmMap()
{
  this->CleanUp();
}
//___________________________________________________________________________
void XSecAlgorithmMap::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void XSecAlgorithmMap::Init(void)
{
  fEventGeneratorList = 0;

  fInitState       = new InitialState;
  fInteractionList = new InteractionList;
}
//___________________________________________________________________________
void XSecAlgorithmMap::CleanUp(void)
{
  delete fInitState;
  delete fInteractionList;

  this->clear();
}
//___________________________________________________________________________
void XSecAlgorithmMap::Copy(const XSecAlgorithmMap & xsmap)
{
  fEventGeneratorList = xsmap.fEventGeneratorList;

  fInitState       -> Copy (*xsmap.fInitState);
  fInteractionList -> Copy (*xsmap.fInteractionList);

  this->clear();

  XSecAlgorithmMap::const_iterator iter;

  for(iter = xsmap.begin(); iter != xsmap.end(); ++iter) {
    string code = iter->first;
    const XSecAlgorithmI * alg = iter->second;

    this->insert(map<string, const XSecAlgorithmI *>::value_type(code,alg));
  }
}
//___________________________________________________________________________
void XSecAlgorithmMap::UseGeneratorList(const EventGeneratorList * list)
{
  fEventGeneratorList = list;
}
//___________________________________________________________________________
void XSecAlgorithmMap::BuildMap(const InitialState & init_state)
{
  LOG("XSecAlgMap", pNOTICE)
                << "Building 'interaction' -> 'xsec algorithm' associations";
  LOG("XSecAlgMap", pNOTICE)
           << "Using all simulated interactions for init-state: "
                                                    << init_state.AsString();
  if(!fEventGeneratorList) {
     LOG("XSecAlgMap", pWARN)
      << "No EventGeneratorList was loaded. Will not build XSecAlgorithmMap";
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
     LOG("XSecAlgMap", pNOTICE)
        << "Querying [" << evgen->Id().Key() << "] for its InteractionList";

     const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
     InteractionList * ilst = ilstgen->CreateInteractionList(init_state);

     // no point to go on if the list is NULL - continue to next iteration
     if(!ilst) continue;

     // append the new InteractionList to the local copy
     fInteractionList->Append(*ilst);

     // cross section algorithm used by this EventGenerator
     const XSecAlgorithmI * xsec_alg = evgen->CrossSectionAlg();

     // loop over all interaction that can be genererated by the current
     // EventGenerator and link all of them to the current XSecAlgorithmI
     for(intliter = ilst->begin(); intliter != ilst->end(); ++intliter)
     {
        // current interaction
        Interaction * interaction = *intliter;
        string code = interaction->AsString();

         // link with the xsec algorithm
         SLOG("XSecAlgMap", pINFO)
           << "\nLinking: " << code
              << "\n     --> with xsec algorithm: " << xsec_alg->Id().Key();
         this->insert(
            map<string, const XSecAlgorithmI *>::value_type(code,xsec_alg));

     } // loop over interactions
     delete ilst;
     ilst = 0;
  } // loop over event generators
}
//___________________________________________________________________________
const XSecAlgorithmI * XSecAlgorithmMap::FindXSecAlgorithm(
                                      const Interaction * interaction) const
{
  if(!interaction) {
    LOG("XSecAlgMap", pWARN) << "Null interaction!!";
    return 0;
  }

  string code = interaction->AsString();

  XSecAlgorithmMap::const_iterator xsec_alg_iter = this->find(code);
  if(xsec_alg_iter == this->end()) {
    LOG("XSecAlgMap", pWARN)
         << "No XSecAlgorithmI was found for interaction: \n" << code;
    return 0;
  }

  const XSecAlgorithmI * xsec_alg = xsec_alg_iter->second;
  return xsec_alg;
}
//___________________________________________________________________________
const InteractionList & XSecAlgorithmMap::GetInteractionList(void) const
{
  return *fInteractionList;
}
//___________________________________________________________________________
void XSecAlgorithmMap::Print(ostream & stream) const
{
  XSecAlgorithmMap::const_iterator iter;

  stream<< "Printing 'interaction' -> 'xsec algorithm' associations" << endl;

  for(iter = this->begin(); iter != this->end(); ++iter) {
    string code = iter->first;
    const XSecAlgorithmI * alg = iter->second;
    if(alg) {
       stream << code << "  -> " << alg->Id().Key() << endl;
    } else {
       stream << code << "  ->  **** NULL XSEC ALGORITHM ****" << endl;
    }
  }
}
//___________________________________________________________________________
XSecAlgorithmMap & XSecAlgorithmMap::operator = (const XSecAlgorithmMap & xs)
{
  this->Copy(xs);
  return (*this);
}
//___________________________________________________________________________


