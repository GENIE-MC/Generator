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

#include "TVector3.h"

namespace genie {

class InitialState;

namespace vmc {

class IRecordList;
class Record;

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
  struct Key
  {
    Key(int _nucl_pdg, int _nu_pdg, bool _iscc)
      : nucl_pdg(_nucl_pdg), nu_pdg(_nu_pdg), iscc(_iscc) {}

    bool operator<(const Key& k) const
    {
      return (std::make_tuple(  nucl_pdg,   nu_pdg,   iscc) <
              std::make_tuple(k.nucl_pdg, k.nu_pdg, k.iscc));
    }

    friend std::ostream& operator<<(std::ostream& os, const Key& k)
    {
      os << k.nu_pdg << " on " << k.nucl_pdg << " " << " via " << (k.iscc ? "CC" : "NC");
      return os;
    }

    int nucl_pdg;
    int nu_pdg;
    bool iscc;
  };

  mutable std::map<Key, const IRecordList*> fRecords;

  const Record* GetRecord(const InitialState& init_state) const;

  /// Return a random (x,y,z) basis with z aligned with the input vector
  std::vector<TVector3> Basis(TVector3 z) const;

  void LoadRecords() const;
};

} // vmc namespace
} // genie namespace

#endif // _EVENT_LIBRARY_INTERFACE_H_
