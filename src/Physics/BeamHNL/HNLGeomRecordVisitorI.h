//____________________________________________________________________________
/*!

\class   genie::hnl::GeomRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the HNL VertexGenerator module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 Costas Andreopoulos <c.andreopoulos \at cern.ch>
	 University of Liverpool

\created February 16th, 2023

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_GEOM_RECORD_VISITOR_I_H_
#define _HNL_GEOM_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

  class GHepRecord;

  namespace hnl {

    class GeomRecordVisitorI: public EventRecordVisitorI {

    public:
      
      virtual ~GeomRecordVisitorI();

      //-- define the GeomRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual void SetGeomFile( std::string geomfile ) const = 0;

    protected:

      GeomRecordVisitorI();
      GeomRecordVisitorI(string name);
      GeomRecordVisitorI(string name, string config);

    };

  } // namespace hnl

} // namespace genie

#endif // #ifndef _HNL_GEOM_RECORD_VISTOR_I_H_
