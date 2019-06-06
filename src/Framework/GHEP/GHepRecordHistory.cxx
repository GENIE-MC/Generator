//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - September 23, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TSystem.h>

#include "Framework/GHEP/GHepRecordHistory.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/PrintUtils.h"

using std::endl;

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const GHepRecordHistory & history)
 {
   history.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
GHepRecordHistory::GHepRecordHistory() :
map<int, GHepRecord*>()
{
  this->ReadFlags();
}
//___________________________________________________________________________
GHepRecordHistory::GHepRecordHistory(const GHepRecordHistory & history) :
map<int, GHepRecord*>()
{
  this->Copy(history);
  this->ReadFlags();
}
//___________________________________________________________________________
GHepRecordHistory::~GHepRecordHistory()
{
  this->PurgeHistory();
}
//___________________________________________________________________________
void GHepRecordHistory::AddSnapshot(int step, GHepRecord * record)
{
// Adds a GHepRecord 'snapshot' at the history buffer

  bool go_on = (fEnabledFull || (fEnabledBootstrapStep && step==-1));
  if(!go_on) return;

  if(!record) {
   LOG("GHEP", pWARN)
    << "Input GHEP record snapshot is null. Is not added at history record";
    return;
  }

  if( this->count(step) == 0 ) {

     LOG("GHEP", pNOTICE)
                     << "Adding GHEP snapshot for processing step: " << step;

     GHepRecord * snapshot = new GHepRecord(*record);
     this->insert( map<int, GHepRecord*>::value_type(step,snapshot));

  } else {
     // If you have already stepped back and reprocessing, then you should
     // have purged the 'recent' history (corresponing to 'after the return
     // processing step')
     LOG("GHEP", pWARN)
      << "GHEP snapshot for processing step: " << step << " already exists!";
  }
}
//___________________________________________________________________________
void GHepRecordHistory::PurgeHistory(void)
{
  LOG("GHEP", pNOTICE) << "Purging GHEP history buffer";

  GHepRecordHistory::iterator history_iter;
  for(history_iter = this->begin();
                              history_iter != this->end(); ++history_iter) {

    int step = history_iter->first;
    LOG("GHEP", pINFO) 
                  << "Deleting GHEP snapshot for processing step: " << step;

    GHepRecord * record = history_iter->second;
    if(record) {
      delete record;
      record = 0;
    }
  }
  this->clear();
}
//___________________________________________________________________________
void GHepRecordHistory::PurgeRecentHistory(int start_step)
{
// Snapshots are added to the history record *after* each processing step
// (marked 0,1,2,...). A special snapshot corresponding to the event record
// before any processing step is added with key = -1.
// Therefore GHepRecordHistory keys should be: -1,0,1,2,3,...

  LOG("GHEP", pNOTICE) 
       << "Purging recent GHEP history buffer (processing step >= " 
                                                      << start_step << ")";

  if(start_step < -1) {
    LOG("GHEP", pWARN) 
               << "Invalid starting step: " << start_step << " - Ignoring";
    return;
  }

  if(start_step == -1) {
    // delete everything
    this->PurgeHistory();
    return;
  }

  GHepRecordHistory::iterator history_iter;
  for(history_iter = this->begin();
                              history_iter != this->end(); ++history_iter) {

    if(history_iter->first >= start_step) { 
       int step = history_iter->first;
       LOG("GHEP", pINFO) 
                  << "Deleting GHEP snapshot for processing step: " << step;
       this->erase(history_iter); 
    }
  }
}
//___________________________________________________________________________
void GHepRecordHistory::Copy(const GHepRecordHistory & history)
{
  this->PurgeHistory();

  GHepRecordHistory::const_iterator history_iter;
  for(history_iter = history.begin();
                           history_iter != history.end(); ++history_iter) {

    unsigned int step   = history_iter->first;
    GHepRecord * record = history_iter->second;

    this->AddSnapshot(step, record);
  }
}
//___________________________________________________________________________
void GHepRecordHistory::Print(ostream & stream) const
{
  stream << "\n ****** Printing GHEP record history"
                              << " [depth: " << this->size() << "]" << endl;

  GHepRecordHistory::const_iterator history_iter;
  for(history_iter = this->begin();
                              history_iter != this->end(); ++history_iter) {

    unsigned int step   = history_iter->first;
    GHepRecord * record = history_iter->second;

    stream << "\n[After processing step = " << step << "] :";

    if(!record) {
      stream
         << "** ERR: No history record available for this processing step!";
    } else {
      stream << *record;
    }
  }
}
//___________________________________________________________________________
void GHepRecordHistory::ReadFlags(void) 
{
  if (gSystem->Getenv("GHEPHISTENABLE")) {

     string envvar = string(gSystem->Getenv("GHEPHISTENABLE"));

     fEnabledFull          = (envvar=="FULL")      ? true:false;
     fEnabledBootstrapStep = (envvar=="BOOTSTRAP") ? true:false;

  } else {
     // set defaults
     fEnabledFull          = false;          
     fEnabledBootstrapStep = true; 
  }

  LOG("GHEP", pINFO) << "GHEP History Flags: ";
  LOG("GHEP", pINFO) << "  - Keep Full History:          " 
                     << utils::print::BoolAsYNString(fEnabledFull);
  LOG("GHEP", pINFO) << "  - Keep Bootstrap Record Only: " 
                     << utils::print::BoolAsYNString(fEnabledBootstrapStep);
}
//___________________________________________________________________________

