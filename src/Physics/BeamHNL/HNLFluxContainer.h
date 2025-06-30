//____________________________________________________________________________
/*!

\class    genie::hnl::FluxContainer

\brief    A GENIE flux container specific for HNL containers.
          Based on the dk2nu flux paradigm and genie::flux::GNuMIFluxPassThroughInfo
	  
	  Also see $GENIE/src/contrib/beamhnl/write_dk2nus.C for an example of 
	  expected flux-input structure

\author   John Plows

\created  Feb 16, 2023

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _HNL_FLUX_CONTAINER_H_
#define _HNL_FLUX_CONTAINER_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/GFluxI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"

#include "Physics/BeamHNL/HNLDecayUtils.h"

class TFile;
class TChain;
class TTree;
class TBranch;

using std::string;
using std::ostream;

namespace genie{
  
  namespace hnl {

    class FluxContainer;
    ostream & operator << (ostream & stream, const FluxContainer & gnmf);

    /// FluxContainer:
    /// =============================
    /// A C-struct that is based on the flux::GNuMIFluxPassThroughInfo
    /// struct, but which is not GNuMI specific. Accepts enough information 
    /// about the HNL fluxes and the base flux to pass all the necessary
    /// information to the hnl::FluxCreator class for flux calculations.
    /// =============================
 
    class FluxContainer: public TObject {
  
    public: 
      FluxContainer();
      virtual ~FluxContainer() {};

      void ResetCopy() const;
      void Print(const Option_t * /* opt */) const;
 
      friend ostream & operator << (ostream & stream, const FluxContainer & gnmf);
 
      // members

      mutable int evtno;                 ///< Event number

      mutable int pdg;                   ///< HNL PDG code
      mutable int parPdg;                ///< parent PDG code
      mutable int lepPdg;                ///< PDG code of lepton produced with HNL on parent decay
      mutable int nuPdg;                 ///< PDG code of SM neutrino that would have been produced

      mutable int prodChan;              ///< Decay mode that produced HNL
      mutable int nuProdChan;            ///< Decay mode that would have produced SM neutrino

      mutable TVector3 startPoint;       ///< parent decay vertex in NEAR coords [m]
      mutable TVector3 targetPoint;      ///< point in detector HNL is forced towards in NEAR coords [m]
      mutable TVector3 startPointUser;   ///< parent decay vertex in USER coords [m]
      mutable TVector3 targetPointUser;  ///< point in detector HNL is forced towards in USER coords [m]
      mutable double delay;              ///< delay HNL would have wrt SMv [ns]

      mutable TVector3 polz;             ///< HNL polarisation vector, in HNL rest frame, in NEAR coords
      
      mutable TLorentzVector p4;         ///< HNL momentum in NEAR coords [GeV/c]
      mutable TLorentzVector parp4;      ///< parent momentum at HNL production in NEAR coords [GeV/c]
      mutable TLorentzVector p4User;     ///< HNL momentum in USER coords [GeV/c]
      mutable TLorentzVector parp4User;  ///< parent momentum at HNL production in USER coords [GeV/c]

      mutable double Ecm;                ///< Parent rest-frame energy of HNL [GeV]
      mutable double nuEcm;              ///< Parent rest-frame energy of equivalent SM neutrino [GeV]
      
      mutable double XYWgt;              ///< geometric acceptance (angular size of detector in parent rest frame)
      mutable double boostCorr;          ///< boost correction wrt parent rest-frame (ELAB = ECM * boostCorr)

      mutable double accCorr;            ///< acceptance correction (collimation effect. SM v == 1)
      mutable double zetaMinus;          ///< minimum angular deviation from parent momentum to reach detector [deg]
      mutable double zetaPlus;           ///< maximum angular deviation from parent momentum to reach detector [deg]

      mutable double acceptance;         ///< full acceptance == XYWgt * boostCorr^2 * accCorr

      mutable double nimpwt;             ///< Weight of parent
      

    }; // class FluxContainer

  } // namespace hnl

} // namespace genie

#endif // #ifndef _HNL_FLUX_CONTAINER_H_
