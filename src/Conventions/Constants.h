//____________________________________________________________________________
/*!

\namespace genie::constants

\brief     Physics constants

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 03, 2004

*/
//____________________________________________________________________________

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <TMath.h>

#include "Conventions/Units.h"

namespace genie {
namespace constants {

//----- Basic constants

static const double kLightSpeed    = 1.;
static const double kPlankConstant = 1.;

//----- pi, e,...

static const double kPi    = 3.1415927;
static const double kPi2   = TMath::Power(kPi,2);
static const double kPi3   = TMath::Power(kPi,3);
static const double kPi4   = TMath::Power(kPi,4);
static const double ke     = 2.7182818;
static const double kSqrte = TMath::Sqrt(ke);

//----- Avogadro number, compton wavelength and such...

static const double kNA    = 6.023e23;             
static const double kLe    = 3.8616E-11 *units::cm;
static const double kLe2   = TMath::Power(kLe,2);          

//----- Coupling constants

static const double kAem   = 1./137.03599976; // dimensionless - EM coupling const
static const double kAem2  = TMath::Power(kAem,2);

static const double kGF    = 1.16639E-5;      // GeV^-2 - Fermi const from b-decay
static const double kGF2   = TMath::Power(kGF,2);

//----- Masses

//      For simplicity, the most commonly used particle masses defined here.
//      In general, however, particle masses in GENIE classes should be obtained
//      through the genie::PDGLibrary as shown below:
//           double mass = PDGLibrary::Instance()->Find(pdg_code)->Mass();
//      For consistency, the values below must match whatever is used in
//      PDGLibrary.

static const double kElectronMass   =  0.0005109989;        // GeV
static const double kMuonMass       =  0.105658357;         // GeV
static const double kTauMass        =  1.77703;             // GeV
static const double kPionMass       =  0.140;               // GeV
static const double kProtonMass     =  0.9382720;           // GeV
static const double kNeutronMass    =  0.9395653;           // GeV
static const double kNucleonMass    =  (kProtonMass+kNeutronMass)/2.;

static const double kElectronMass2  =  TMath::Power(kElectronMass,2); // GeV^2
static const double kMuonMass2      =  TMath::Power(kMuonMass,2);     // GeV^2
static const double kTauMass2       =  TMath::Power(kTauMass,2);      // GeV^2
static const double kPionMass2      =  TMath::Power(kPionMass,2);     // GeV^2
static const double kProtonMass2    =  TMath::Power(kProtonMass,2);   // GeV^2
static const double kNeutronMass2   =  TMath::Power(kNeutronMass,2);  // GeV^2
static const double kNucleonMass2   =  TMath::Power(kNucleonMass,2);  // GeV^2

static const double kMw             =  80.14;                // GeV - W boson mass
static const double kMz             =  91.19;                // GeV - Z boson mass
static const double kMw2            =  TMath::Power(kMw,2);  // GeV^2
static const double kMz2            =  TMath::Power(kMz,2);  // GeV^2

//----- sqrts frequently encountered in Helicity Amplitude calculations

static const double kSqrt2  = TMath::Sqrt(  2.0 );
static const double kSqrt3  = TMath::Sqrt(  3.0 );
static const double kSqrt4  = 2.0;
static const double kSqrt5  = TMath::Sqrt(  5.0 );
static const double kSqrt6  = TMath::Sqrt(  6.0 );
static const double kSqrt7  = TMath::Sqrt(  7.0 );
static const double kSqrt8  = TMath::Sqrt(  8.0 );
static const double kSqrt9  = TMath::Sqrt(  9.0 );
static const double kSqrt10 = TMath::Sqrt( 10.0 );
static const double kSqrt15 = TMath::Sqrt( 15.0 );
static const double kSqrt18 = TMath::Sqrt( 18.0 );
static const double kSqrt24 = TMath::Sqrt( 24.0 );
static const double kSqrt27 = TMath::Sqrt( 27.0 );
static const double kSqrt30 = TMath::Sqrt( 30.0 );
static const double kSqrt35 = TMath::Sqrt( 35.0 );

//----- Misc constants for empirical formulas

// Ro in nuclear radius formula R=Ro*A^(1/3), in GeV^-1
static const double kNucRo = 1.2E-15*units::m;
// Nuclear density (in nuclear core), in GeV^4
static const double kNucDensity = 2.3E+17 *units::kg/units::m3;

} // namespace constants
} // namespace genie

#endif // _CONSTANTS_H_


