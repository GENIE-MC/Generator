//____________________________________________________________________________
/*!

\class   genie::hnl::FluxRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the HNL FluxCreator module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	 University of Liverpool & STFC Rutherford Appleton Laboratory

\created January 22nd, 2023

\cpright Copyright (c) 2003-2022, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _GNUMI_FLUX_RECORD_VISITOR_I_H_
#define _GNUMI_FLUX_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

  class GHepRecord;

  namespace flux {
    class GNuMIFluxPassThroughInfo;
  }

  namespace hnl {

    class FluxRecordVisitorI: public EventRecordVisitorI {

    public:

      virtual ~FluxRecordVisitorI();

      //-- define the FluxRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual flux::GNuMIFluxPassThroughInfo * RetrieveGNuMIFluxPassThroughInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxInfo() const = 0;
      virtual flux::GNuMIFluxPassThroughInfo RetrieveFluxBase() const = 0;

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

#endif // #ifndef _GNUMI_FLUX_RECORD_VISITOR_I_H_
