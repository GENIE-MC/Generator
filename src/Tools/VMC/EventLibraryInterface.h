//____________________________________________________________________________
/*!

\class    genie::vmc::EventLibraryInterface

\brief    Reads pre-generated events produced by an external generator.
          On an event-by-event basis, it can accept GENIE input specifying
          the neutrino and target IDs and neutrino energy and, therefore,
          it can re-use the upstream GENIE flux and geometry tools.

\author   

\created  February 28, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _EVENT_LIBRARY_INTERFACE_H_
#define _EVENT_LIBRARY_INTERFACE_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {
namespace vmc {

class EventLibraryInterface: public EventRecordVisitorI {

public:
  EventLibraryInterface();
  EventLibraryInterface(string config);
 ~EventLibraryInterface();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

   void LoadConfig            (void);
   void AddInitialState       (GHepRecord * event) const;
   void GenerateDecayProducts (GHepRecord * event) const;
};

} // vmc namespace
} // genie namespace

#endif // _EVENT_LIBRARY_INTERFACE_H_
