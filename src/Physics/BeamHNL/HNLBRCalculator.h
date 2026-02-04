//----------------------------------------------------------------------------
/*!

  An Algorithm that manages the calculations of production and decay
  branching ratios for HNL.

\class      genie::hnl::BRCalculator

\brief      Manages HNL BR (prod and decay)

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 3rd, 2021

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _HNL_BRFUNCTIONS_H_
#define _HNL_BRFUNCTIONS_H_

// -- C++ includes
#include <iterator>
#include <sstream>

// -- ROOT includes
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

// -- GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

#include "Physics/BeamHNL/HNLChannelCalculatorI.h"

#include "Physics/BeamHNL/FormFactorTables.h"
#include "Physics/BeamHNL/HNLDecayUtils.h"
#include "Physics/BeamHNL/HNLKinUtils.h"
#include "Physics/BeamHNL/HNLProductionMode.h"

namespace genie {

  class PDGLibrary;
  
  namespace hnl {

    class BRCalculator : public ChannelCalculatorI {

    public:

      BRCalculator();
      BRCalculator(string name);
      BRCalculator(string name, string config);
      ~BRCalculator();

      // -- implement the ChannelCalculatorI interface

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

      // return the kinematic scaling for a production channel
      double KinematicScaling( genie::hnl::HNLProd_t hnlprod ) const;
      
      // return the integrated decay width for a decay channel
      double DecayWidth( genie::hnl::HNLDecayMode_t hnldm ) const;

    private:

      void LoadConfig(void);

      //============================================
      // total decay widths, parents to HNL
      //============================================
      
      double KScale_Global( genie::hnl::HNLProd_t hnldm, const double M ) const;
      
      // P --> N + \ell_{\alpha} 
      double DWidth_PseudoscalarToLepton( const double mP, const double M, const double Ua42, const double ma ) const;
      double KScale_PseudoscalarToLepton( const double mP, const double M, const double ma ) const; // uses formulae (2.12)-(2.14) from Shrock, Phys. Rev. D 24 (1981) 5
	
      // P --> N + \ell_{\alpha} + \pi^0
      double DWidth_PseudoscalarToPiLepton( const double mP, const double M, const double Ua42, const double ma ) const;
      double KScale_PseudoscalarToPiLepton( const double mP, const double M, const double ma ) const;

      // mu --> N + \nu_{\alpha} + e
      double DWidth_MuonToNuAndElectron( const double M, const double Ue42, const double Umu42, const double Ut42 ) const;
      double KScale_MuonToNuAndElectron( const double M ) const;

      //============================================
      // total decay widths for HNL channels
      //============================================

      double DWidth_Global( genie::hnl::HNLDecayMode_t hnldm, const double M ) const;
	
      // N --> pi0 + nu_\alpha
      double DWidth_PiZeroAndNu( const double M, const double Ue42, const double Umu42, const double Ut42 ) const;
	
      // N --> pi^{\pm} + \ell_{\alpha}
      double DWidth_PiAndLepton( const double M, const double Ua42, const double ma ) const;

      // invisible
      double DWidth_Invisible( const double M, const double Ue42, const double Umu42, const double Ut42 ) const;

      // N --> nu + \ell_{\beta} \ell_{\beta} , \beta = e or mu
      double DWidth_SameLepton( const double M, const double Ue42, const double Umu42, const double Ut42, const double mb, bool bIsMu ) const;

      // N --> nu + \ell_{\alpha}^{-} + \ell_{\beta}^{+}, \alpha != \beta, \alpha = e or mu
      // alpha is the "correct" lepton number sign ( Q < 0 for N, Q > 0 for Nbar)
      double DWidth_DiffLepton( const double M, const double Ua42, const double Ub42, const int IsMajorana ) const;

      // N --> pi^\pm + pi^0 + \ell_{\alpha}^\mp, \alpha = e or mu
      double DWidth_PiPi0Ell( const double M, const double ml,
			      const double Ue42, const double Umu42, const double Ut42,
			      const bool isElectron = false ) const;

      // N --> pi^0 + pi^0 + \nu (all flavours summed up)
      double DWidth_Pi0Pi0Nu( const double M,
			      const double Ue42, const double Umu42, const double Ut42 ) const;

      static double PiPi0EllForm( double *x, double *par );
      static double Pi0Pi0NuForm( double *x, double *par );

      // kinematic functions
      double GetFormfactorF1( double x ) const;
      double GetFormfactorF2( double x ) const;
     
      bool fIsConfigLoaded = false;
      
      // physical constants
      double wAng, s2w;
      const double GF  = genie::constants::kGF; // GeV^{-2}
      const double GF2 = GF*GF;
      const double pi  = genie::constants::kPi;
      double fpi, fpi2;
      double BR_C1, BR_C2;
		
      double mPi0, mPi, mMu, mK, mK0, mE;
      double fMass;
      std::vector< double > fCouplings; double fUe42, fUm42, fUt42;
      bool fMajorana;

      double Vud, Vud2;

      // PMNS matrix elements
      double Ue1, Ue2, Ue3, Um1, Um2, Um3, Ut1, Ut2, Ut3;

      // scaling factors (3-body decays) from Ballett et al, JHEP 2020 (2019) 3 Fig. 2
      // this figure was digitised so write as map between x = HNL mass [GeV] and y = scaling
      // Digitisation done using WebPlotDigitizer (https://apps.automeris.io/wpd ; https://github.com/ankitrohatgi/WebPlotDigitizer)
      std::map< double, double > kscale_K3mu, kscale_K3e, kscale_mu3e;
    
    }; // class BRCalculator

} // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_BRFUNCTIONS_H_
