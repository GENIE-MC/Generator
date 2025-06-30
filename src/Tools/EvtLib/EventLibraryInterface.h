//____________________________________________________________________________
/*!

\class    genie::evtlib::EventLibraryInterface

\brief    Reads pre-generated events produced by an external generator.
          On an event-by-event basis, it can accept GENIE input specifying
          the neutrino and target IDs and neutrino energy and, therefore,
          it can re-use the upstream GENIE flux and geometry tools.

\author   

\created  February 28, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _EVENT_LIBRARY_INTERFACE_H_
#define _EVENT_LIBRARY_INTERFACE_H_

#include "Tools/EvtLib/Key.h"
class IEvtLibRecordList;

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Kinematics.h"

#include "TVector3.h"
#include "TFile.h"

#include <map>

namespace genie {

class Interaction;

namespace evtlib {

class IEvtLibRecordList;
class EvtLibRecord;

class EventLibraryInterface: public EventRecordVisitorI {

public:
  EventLibraryInterface();
  EventLibraryInterface(string config);
 ~EventLibraryInterface();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const override ;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config) override ;
  void Configure(string config) override ;

 protected:

  const EvtLibRecord* GetRecord(const Interaction* interaction) const;

  void LoadRecords();
  void Cleanup();

  void FillKinematics( const GHepRecord &,
                       genie::Kinematics & kine,
                       int primary_lep_id ) const ;

private:

  std::map<Key, const IEvtLibRecordList*> fRecords;
  TFile* fRecordFile;
};

} // evtlib namespace
} // genie namespace

#endif // _EVENT_LIBRARY_INTERFACE_H_
