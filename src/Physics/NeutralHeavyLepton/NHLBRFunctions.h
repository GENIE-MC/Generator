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
// TODO: add in decays to ~pi pi0 ell~ and pi0 pi0 nuell
//       Figure out how to read Ue4, Umu4, Ut4 from config
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
#include "Physics/NeutralHeavyLepton/NHLKinUtils.h"

namespace genie {
namespace NHL {

    namespace NHLSelector {

	// physical constants, PDG 2021
	static const double s2w = 0.22767; // Weinberg mixing angle consistent with genie kMw, kMz
	static const double GF  = genie::constants::kGF; // GeV^{-2}
	static const double GF2 = GF*GF;
	static const double pi  = genie::constants::kPi;
	static const double fpi = 0.130 * genie::units::GeV; // GeV
	static const double fpi2 = fpi*fpi;
    
	static const double mPi0 = genie::constants::kPi0Mass;
	static const double mPi  = genie::constants::kPionMass;
	static const double mMu  = genie::constants::kMuonMass;

	static const double Vud  = 0.9737; // PDG 2021
	static const double Vud2 = Vud * Vud;

	// PMNS matrix elements, assumed real
	// using CV values from v5.1 of www.nu-fit.org - accessed Jan 25th 2022
	// using SK-atm data but ignoring correlations!
	static const double Ue1 = 1./2. * ( 0.801 + 0.845 );
	static const double Ue2 = 1./2. * ( 0.513 + 0.579 );
	static const double Ue3 = 1./2. * ( 0.144 + 0.156 );
	static const double Um1 = 1./2. * ( 0.244 + 0.499 );
	static const double Um2 = 1./2. * ( 0.505 + 0.693 );
	static const double Um3 = 1./2. * ( 0.631 + 0.763 );
	static const double Ut1 = 1./2. * ( 0.272 + 0.518 );
	static const double Ut2 = 1./2. * ( 0.471 + 0.669 );
	static const double Ut3 = 1./2. * ( 0.623 + 0.761 );
    
	static const double BR_C1 = 1./4. * ( 1. - 4. * s2w + 8. * s2w * s2w );
	static const double BR_C2 = 1./2. * ( -s2w + 2. * s2w * s2w );

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
