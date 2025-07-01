//____________________________________________________________________________
/*!

\class    genie::BSKLNBaseRESPXSec2014

\brief    Base class for the Berger-Sehgal and the Kuzmin, Lyubushkin, Naumov
          resonance models, implemented as modifications to the Rein-Sehgal model.

\ref      Berger, Sehgal Phys. Rev. D76, 113004 (2007) \n
          Kuzmin, Lyubushkin, Naumov Mod. Phys. Lett. A19 (2004) 2815 \n

          D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances

          and Single Pion Production, Ann.Phys.133, 79 (1981) \n

          
          K. M. Graczyk* and J. T. Sobczyk,PHYSICAL REVIEW D 77, 053001 (2008) for vector and axial current calculation \n

          Modifications based on a MiniBooNE tune courtesy of J. Nowak, S.Dytman

\author   Steve Dytman
          University of Pittsburgh

          Jarek Nowak
          University of Lancaster

          Gabe Perdue
          Fermilab
          
          Igor Kakorin <kakorin@inr.ru>
          Joint Institute for Nuclear Research 
          

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Sep 15, 2015

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _BSKLN_BASE_RES_PXSEC_2014_H_
#define _BSKLN_BASE_RES_PXSEC_2014_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/Resonance/XSection/FKR.h"

namespace genie {

  class RSHelicityAmplModelI;
  class XSecIntegratorI;

  class BSKLNBaseRESPXSec2014: public XSecAlgorithmI {

    public:
      virtual ~BSKLNBaseRESPXSec2014();

      // implement the XSecAlgorithmI interface
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    protected:

      BSKLNBaseRESPXSec2014(string name);
      BSKLNBaseRESPXSec2014(string name, string config);

      void LoadConfig (void);

      mutable FKR fFKR;

      const RSHelicityAmplModelI * fHAmplModelCC;
      const RSHelicityAmplModelI * fHAmplModelNCp;
      const RSHelicityAmplModelI * fHAmplModelNCn;
      const RSHelicityAmplModelI * fHAmplModelEMp;
      const RSHelicityAmplModelI * fHAmplModelEMn;

      double fFermiConstant2 ;
      double fFineStructure2 ;

      // configuration data
      bool     fWghtBW;            ///< weight with resonance breit-wigner?
      bool     fNormBW;            ///< normalize resonance breit-wigner to 1?
      double   fZeta;              ///< FKR parameter Zeta
      double   fOmega;             ///< FKR parameter Omega
      double   fCa50;              ///< CA5_0
      double   fMa2;               ///< (axial mass)^2
      double   fMv2;               ///< (vector mass)^2
      double   fVud2;              ///< |Vud|^2(square of magnitude ud-element of CKM-matrix)
      bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
      double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
      double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
      double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
      double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
      string fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
      bool fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
      bool fUsePauliBlocking;      ///< account for Pauli blocking?

      double   fXSecScaleCC;       ///< external CC xsec scaling factor
      double   fXSecScaleNC;       ///< external NC xsec scaling factor
      double   fXSecScaleEM;       ///< external EM xsec scaling factor

      bool fKLN;
      bool fBRS;

      // Use tuned axial, vector form factors by MiniBooNE
      // Tuned to ANL BNL data
      bool fGAMiniBooNE;
      bool fGVMiniBooNE;

      // GV calculation coeff
      double fCv3;
      double fCv4;
      double fCv51;
      double fCv52;

      // Sarita-Schwinger prenscription parameters from PhysRevD.77.053001
      bool   fGASaritaSchwinger ; 
      bool   fGVSaritaSchwinger ; 
      double fcII ; 
      double fMb2 ; 

      const XSecIntegratorI * fXSecIntegrator;
  };

}       // genie namespace

#endif  // _BSKLN_BASE_RES_PXSEC_2014_H_
