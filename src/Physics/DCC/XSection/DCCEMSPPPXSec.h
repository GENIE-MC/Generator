//____________________________________________________________________________
/*!

\class    genie::DCCEMSPPPXSec

\brief    Class calculate differental cross-section d^4(sig)/d(Q2)d(W)d(cost)d(phi),
          for specific W, Q2, neutrino energy(in lab frame) & pion angles in the Adler frame, where \n
          Q2          : Sqaured 4-momentum transfer, Q2 = -k*k         \n                        
          W           : Invariant mass                                 \n                      
          cost        : Cosine of pion polar angle in Adler frame      \n                      
          phi         : Pion azimuthal angle in Adler frame            \n

for the following channels:      

          1       l + p -> l + p + pi0
          2       l + p -> l + n + pi+
          3       l + n -> l + n + pi0
          4       l + n -> l + p + pi-
                                                                                           
                                                          
\ref      
          1. T. Sato , T.-S. H. Lee, Phys. Rev. C 54, 2660(1996)
          2. A. Matsuyama , T.-S. H. Lee, T. Sato, Phys. Rept. 439, 193(2007)
          3. https://www.phy.anl.gov/theory/research/anl-osaka-pwa/ (and other references in it)
    

\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n

\created  Sep 06, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DCC_EM_SPP_PXSEC_H_
#define _DCC_EM_SPP_PXSEC_H_

#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <type_traits>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResList.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/Resonance/XSection/FKR.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"

namespace genie {

  
  class XSecIntegratorI;

  class DCCEMSPPPXSec: public XSecAlgorithmI {

    
    public:
      DCCEMSPPPXSec();
      DCCEMSPPPXSec(string config);
      virtual ~DCCEMSPPPXSec();

      // implement the XSecAlgorithmI interface 
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;
      bool   ValidKinematics(const Interaction * interaction) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
            
      void LoadConfig (void);
      
      // configuration data
      double   fFermiConstant ;
      double   fCA50;              ///< FKR parameter Zeta
      double   fOmega;             ///< FKR parameter Omega
      double   fMa2;               ///< (axial mass)^2
      double   fMv2;               ///< (vector mass)^2
      double   fSin2Wein;          ///< sin^2(Weingberg angle)
      double   fVud;               ///< |Vud| (magnitude ud-element of CKM-matrix)
      double   fXSecScaleCC;       ///< external CC xsec scaling factor
      double   fXSecScaleNC;       ///< external NC xsec scaling factor
      const ELFormFactorsModelI  * fElFFModel;          ///< Elastic form facors model for background contribution
      const QELFormFactorsModelI * fFormFactorsModel;   ///< Quasielastic form facors model for background contribution
      const QELFormFactorsModelI * fEMFormFactorsModel; ///< Electromagnetic form factors model for background contribution
      
      string fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
      bool fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
      bool fUsePauliBlocking;      ///< account for Pauli blocking?

      mutable QELFormFactors  fFormFactors;      ///<  Quasielastic form facors for background contribution
      mutable QELFormFactors  fEMFormFactors;    ///<  Electromagnetic form facors for background contribution
      double  f_pi;                              ///<  Constant for pion-nucleon interaction
      double  FA0;                               ///<  Axial coupling (value of axial form factor at Q2=0)
      double  Frho0;                             ///<  Value of form fator F_rho at t=0 
      /// Parameters for vector virtual form factor
      /// for background contribution, which equal to:  
      /// 1,                                              W<VWmin   
      /// V3*W^3+V2*W^2+V1*W+V0                     VWmin<W<VWmax
      /// 0                                               W>VWmax
      double fBkgVWmin;
      double fBkgVWmax; 
      double fBkgV3;  
      double fBkgV2;   
      double fBkgV1;   
      double fBkgV0;   
      double fRho770Mass;                        ///< Mass of rho(770) meson 
      
      const XSecIntegratorI * fXSecIntegrator;
      
      BaryonResList  fResList;
                 
  };
  
  
}       // genie namespace

#endif  // _DCC_EM_SPP_PXSEC_H_
