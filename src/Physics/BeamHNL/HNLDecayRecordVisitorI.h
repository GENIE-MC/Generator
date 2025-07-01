//____________________________________________________________________________
/*!

\class   genie::hnl::DecayRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the HNL Decayer module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 Costas Andreopoulos <c.andreopoulos \at cern.ch>
	 University of Liverpool

\created January 23rd, 2023

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_DECAY_RECORD_VISITOR_I_H_
#define _HNL_DECAY_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

  class GHepRecord;

  namespace hnl {

    class DecayRecordVisitorI: public EventRecordVisitorI {

    public:

      virtual ~DecayRecordVisitorI();

      //-- define the DecayRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual double GetHNLLifetime() const = 0;
      virtual double GetHNLMass() const = 0;
      virtual std::vector< double > GetHNLCouplings() const = 0;
      virtual bool IsHNLMajorana() const = 0;

      virtual std::string GetHNLInterestingChannels() const = 0;

      //-- additional methods for particle-gun

      virtual std::vector< double > GetPGunOrigin() const = 0;
      virtual std::vector< double > GetPGunDOrigin() const = 0;
      
      virtual double GetPGunEnergy() const = 0;
      virtual std::vector< double > GetPGunDirection() const = 0;
      virtual std::vector< double > GetPGunDeviation() const = 0;

    protected:

      DecayRecordVisitorI();
      DecayRecordVisitorI(string name);
      DecayRecordVisitorI(string name, string config);

    };

  } // namespace hnl

} // namespace genie

#endif // #ifndef _HNL_DECAY_RECORD_VISITOR_I_H_
