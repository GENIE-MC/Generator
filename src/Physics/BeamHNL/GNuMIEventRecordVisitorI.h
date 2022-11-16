//____________________________________________________________________________
/*!

\class   genie::hnl::GNuMIEventRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interface
         for flux::GNuMIFluxPassThroughInfo objects.
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

    class GNuMIEventRecordVisitorI: public EventRecordVisitorI {

    public:

      virtual ~GNuMIEventRecordVisitorI();

      //-- define the GNuMIEventRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual flux::GNuMIFluxPassThroughInfo * RetrieveGNuMIFluxPassThroughInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxBase() const = 0;

    protected:

      GNuMIEventRecordVisitorI();
      GNuMIEventRecordVisitorI(string name);
      GNuMIEventRecordVisitorI(string name, string config);

    };

  } // namespace hnl

} // namespace genie

#endif // #ifndef _GNUMI_EVENT_RECORD_VISITOR_I_H_
