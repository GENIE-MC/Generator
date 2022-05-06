//----------------------------------------------------------------------------
/*!

  Definitions of the auxiliary functions needed
  to calculate the decay widths of NHL to various channels

\namespace  genie::NHL::NHLSelector

\brief      Decay widths of NHL

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 3rd, 2021

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------
// TODO: recheck decays to ~pi pi0 ell~ and pi0 pi0 nuell
//----------------------------------------------------------------------------

#ifndef _NHL_BRFUNCTIONS_H_
#define _NHL_BRFUNCTIONS_H_

// -- C++ includes
#include <TMath.h>

// -- ROOT includes
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

// -- GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"

#include "Physics/NeutralHeavyLepton/ColomaTables.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
#include "Physics/NeutralHeavyLepton/NHLKinUtils.h"

namespace genie {
  
  namespace NHL {

    namespace NHLSelector {

        // are parameters initialised?
        static __attribute__((unused)) bool fParamsInitialised = false;

	// physical constants

	// RETHERE why does this cause a segfault?
        //static const double wAng = genie::utils::nhl::GetCfgDouble( "Param", "WeakInt", "WeinbergAngle" );
        static __attribute__((unused)) double wAng, s2w;
	static __attribute__((unused)) const double GF  = genie::constants::kGF; // GeV^{-2}
	static __attribute__((unused)) const double GF2 = GF*GF;
	static __attribute__((unused)) const double pi  = genie::constants::kPi;
	static __attribute__((unused)) double fpi, fpi2;
	static __attribute__((unused)) double BR_C1, BR_C2;
		
	static __attribute__((unused)) double mPi0, mPi, mMu, mK, mK0, mE;

	static __attribute__((unused)) double Vud, Vud2;

	// PMNS matrix elements
	static __attribute__((unused)) double Ue1, Ue2, Ue3, Um1, Um2, Um3, Ut1, Ut2, Ut3;

	// initialise these parameters from config files!
	void InitParameters();

	// declaration of PiPi0EllForm as extern?
	extern double PiPi0EllForm( double *x, double *par );
	extern double Pi0Pi0NuForm( double *x, double *par );

	// kinematic functions
	double GetColomaF1( double x );
	double GetColomaF2( double x );

	//============================================
	// total decay widths for NHL channels
	//============================================
	
	// N --> pi0 + nu_\alpha
	double DWidth_PiZeroAndNu( const double M, const double Ue42, const double Umu42, const double Ut42 );
	
	// N --> pi^{\pm} + \ell_{\alpha}
	double DWidth_PiAndLepton( const double M, const double Ua42, const double ma );

	// invisible
	double DWidth_Invisible( const double M, const double Ue42, const double Umu42, const double Ut42 );

	// N --> nu + \ell_{\beta} \ell_{\beta} , \beta = e or mu
	double DWidth_SameLepton( const double M, const double Ue42, const double Umu42, const double Ut42, const double mb, bool bIsMu );

	// N --> nu + \ell_{\alpha}^{-} + \ell_{\beta}^{+}, \alpha != \beta, \alpha = e or mu
	// alpha is the "correct" lepton number sign ( Q < 0 for N, Q > 0 for Nbar)
	double DWidth_DiffLepton( const double M, const double Ua42, const double Ub42, const int IsMajorana );

	// N --> pi^\pm + pi^0 + \ell_{\alpha}^\mp, \alpha = e or mu
	double DWidth_PiPi0Ell( const double M, const double ml,
				const double Ue42, const double Umu42, const double Ut42,
				const bool isElectron = false );

	// N --> pi^0 + pi^0 + \nu (all flavours summed up)
	double DWidth_Pi0Pi0Nu( const double M,
				const double Ue42, const double Umu42, const double Ut42 );

	//============================================
	// differential decay widths for NHL channels
	// want to *sample* kinematics from them so not double rtype
	//============================================

	// N --> pi^{\pm} + \ell_{\alpha}
	// d\Gamma / dcos(\theta), \theta = angle between polVec and pl
	void Diff1Width_PiAndLepton_CosTheta( const double M, const double Ua42,
					      const double ml,
					      double &thePreFac,
					      double &theCnstPart,
					      double &thePropPart );
    
    } // namespace NHLSelector

} // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_BRFUNCTIONS_H_
