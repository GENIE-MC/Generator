//____________________________________________________________________________
/*!

\class   genie::hnl::HNLEventRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the HNL modules.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	 University of Liverpool & STFC Rutherford Appleton Laboratory

\created November 15, 2022

\cpright Copyright (c) 2003-2022, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _GNUMI_EVENT_RECORD_VISITOR_I_H_
#define _GNUMI_EVENT_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

  class GHepRecord;

  namespace flux {
    class GNuMIFluxPassThroughInfo;
  }

  namespace hnl {

    class HNLEventRecordVisitorI: public EventRecordVisitorI {

    public:

      virtual ~HNLEventRecordVisitorI();

      //-- define the HNLEventRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual flux::GNuMIFluxPassThroughInfo * RetrieveGNuMIFluxPassThroughInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxBase() const = 0;

      virtual std::vector< double > GetB2UTranslation() const = 0;
      virtual std::vector< double > GetB2URotation() const = 0;
      virtual std::vector< double > GetDetOffset() const = 0;
      virtual std::vector< double > GetDetRotation() const = 0;
      
      virtual void SetInputFluxPath( std::string finpath ) const = 0;
      virtual int GetNFluxEntries() const = 0;
      virtual void SetGeomFile( std::string geomfile ) const = 0;

    protected:

      HNLEventRecordVisitorI();
      HNLEventRecordVisitorI(string name);
      HNLEventRecordVisitorI(string name, string config);

    };

  } // namespace hnl

} // namespace genie

#endif // #ifndef _GNUMI_EVENT_RECORD_VISITOR_I_H_
