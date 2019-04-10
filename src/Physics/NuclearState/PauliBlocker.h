//____________________________________________________________________________
/*!

\class    genie::PauliBlocker

\brief    Examines whether the generated event should be Pauli blocked.
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

          Other code improvements, additions, fixes were installed by
          Joe Johnston (Univ of Pittsburgh)

\created  October 08, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PAULI_BLOCKER_H_
#define _PAULI_BLOCKER_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Target.h"

namespace genie {

class FermiMomentumTable;

class PauliBlocker : public EventRecordVisitorI {

public :
  PauliBlocker();
  PauliBlocker(string config);
 ~PauliBlocker();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  /// Get the Fermi momentum needed to check Pauli blocking
  double GetFermiMomentum(const Target& tgt, int pdg_Nf,
    double radius = 0.0) const;

private:
   void LoadModelType(void);

   bool fLFG;
   const FermiMomentumTable * fKFTable;
   string fKFTableName;
};

}      // genie namespace

#endif // _PAULI_BLOCKER_H_
