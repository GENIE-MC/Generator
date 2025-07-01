//____________________________________________________________________________
/*!

\class   genie::hnl::FluxRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the HNL FluxCreator module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 Costas Andreopoulos <c.andreopoulos \at cern.ch>
	 University of Liverpool

\created January 22nd, 2023

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_FLUX_RECORD_VISITOR_I_H_
#define _HNL_FLUX_RECORD_VISITOR_I_H_

//#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/BeamHNL/HNLGeomRecordVisitorI.h"
#include "Physics/BeamHNL/HNLFluxContainer.h"

namespace genie {

  class GHepRecord;

  namespace hnl {

    class FluxContainer;

    class FluxRecordVisitorI: public GeomRecordVisitorI {

    public:

      virtual ~FluxRecordVisitorI();

      //-- define the FluxRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual FluxContainer RetrieveFluxInfo() const = 0;

      virtual std::vector< double > GetB2UTranslation() const = 0;
      virtual std::vector< double > GetB2URotation() const = 0;
      virtual std::vector< double > GetDetOffset() const = 0;
      virtual std::vector< double > GetDetRotation() const = 0;
      
      virtual void SetInputFluxPath( std::string finpath ) const = 0;
      virtual void SetGeomFile( std::string geomfile ) const = 0;
      virtual int GetNFluxEntries() const = 0;
      virtual void SetFirstFluxEntry( int iFirst ) const = 0;

    protected:

      FluxRecordVisitorI();
      FluxRecordVisitorI(string name);
      FluxRecordVisitorI(string name, string config);

    };

  } // namespace hnl

} // namespace genie

#endif // #ifndef _HNL_FLUX_RECORD_VISITOR_I_H_
