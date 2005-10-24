//____________________________________________________________________________
/*!

\class   genie::NucBindEnergyAggregator

\brief   A nuclear binding energy 'collector' which visits the event record,
         finds nucleons originating from within a nuclei and subtracts the
         binding energy they had in the nucleus.

         To record this action in the event record a hypothetical BINDINO is
         added to the event record.

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 19, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_BINDING_ENERGY_AGGREGATOR_H_
#define _NUCLEAR_BINDING_ENERGY_AGGREGATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class GHepParticle;

class NucBindEnergyAggregator : public EventRecordVisitorI {

public :

  NucBindEnergyAggregator();
  NucBindEnergyAggregator(string config);
  ~NucBindEnergyAggregator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  GHepParticle * FindMotherNucleus(int ipos, GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _NUCLEAR_BINDING_ENERGY_AGGREGATOR_H_
