//____________________________________________________________________________
/*!

\namespace genie::constants

\brief     Constants

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 03, 2004

*/
//____________________________________________________________________________

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <TMath.h>

#include "Conventions/Units.h"

using namespace genie::units;

namespace genie {
namespace constants {

//----- Basic constants

static const double kLightSpeed    = 1.;
static const double kPlankConstant = 1.;
static const double kPi            = 3.1415927;
static const double kPi_2          = TMath::Power(kPi,2);
static const double kPi_3          = TMath::Power(kPi,3);
static const double kPi_4          = TMath::Power(kPi,4);
static const double ke             = 2.7182818;
static const double kNA            = 6.023e23;    // Avogadro's number

//----- Coupling constants

static const double kAem   = 1./137.03599976; // dimensionless - EM coupling const
static const double kAem_2 = TMath::Power(kAem,2);

static const double kGF    = 1.16639E-5;      // GeV^-2 - Fermi const from b-decay
static const double kGF_2  = TMath::Power(kGF,2);

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

static const double kElectronMass_2 =  TMath::Power(kElectronMass,2); // GeV^2
static const double kMuonMass_2     =  TMath::Power(kMuonMass,2);     // GeV^2
static const double kTauMass_2      =  TMath::Power(kTauMass,2);      // GeV^2
static const double kPionMass_2     =  TMath::Power(kPionMass,2);     // GeV^2
static const double kProtonMass_2   =  TMath::Power(kProtonMass,2);   // GeV^2
static const double kNeutronMass_2  =  TMath::Power(kNeutronMass,2);  // GeV^2
static const double kNucleonMass_2  =  TMath::Power(kNucleonMass,2);  // GeV^2

static const double kMw             =  80.14;                // GeV - W boson mass
static const double kMz             =  91.19;                // GeV - Z boson mass
static const double kMw_2           =  TMath::Power(kMw,2);  // GeV^2
static const double kMz_2           =  TMath::Power(kMz,2);  // GeV^2

//----- Cabbibo & Weinberg angles angles

static const double kCos8c    =  0.97437;                 // cos   (Cabbibo angle)
static const double kCos8c_2  =  TMath::Power(kCos8c,2);  // cos^2 (Cabbibo angle)
static const double kSin8c_2  =  1-kCos8c_2;              // sin^2 (Cabbibo angle)
static const double kSin8c    =  TMath::Sqrt(kSin8c_2);   // sin   (Cabbibo angle)

static const double kSin8w_2  =  0.23117;          // sin^2 (Weinberg angle)
static const double kCos8w_2  =  1-kSin8w_2;       // cos^2 (Weinberg angle)
static const double kSin8w    =  TMath::Sqrt(kSin8w_2);   // sin   (Weinberg angle)
static const double kCos8w    =  TMath::Sqrt(kCos8w_2);   // cos   (Weinberg angle)

//----- CKM elements

static const double kVud    =  0.9738;      // +/-0.0005
static const double kVus    =  0.2200;      // +/-0.0026
static const double kVcd    =  0.224;       // +/-0.012
static const double kVcs    =  0.996;       // +/-0.013
static const double kVud_2  =  TMath::Power(kVud,2);
static const double kVus_2  =  TMath::Power(kVus,2);
static const double kVcd_2  =  TMath::Power(kVcd,2);
static const double kVcs_2  =  TMath::Power(kVcs,2);

//----- Proton & Neutron Anomalous Magnetic Moments

static const double kMuP =  2.7928473; // proton : anomalous magnetic moment
static const double kMuN = -1.913042;  // neutron: anomalous magnetic moment

//----- Axial-vector & Vector Masses & Axial Form Factor at Q^2 = 0

//Elastic:
static const double kElMv         =  0.840;  // Elastic 'vector mass' in GeV
static const double kElMa         =  1.032;  // Elastic 'axial mass' in GeV
static const double kElEta        =  0.;     // Elastic factor eta
static const double kElGa0        =  1.26;   // Elastic factor ga(0)
static const double kElMv2        =  TMath::Power(kElMv,2);
static const double kElMa2        =  TMath::Power(kElMa,2);

//Quasi-Elastic:
static const double kQelMv        =  0.840;  // Quasi-Elastic 'vector mass' in GeV
static const double kQelMa        =  1.032;  // Quasi-Elastic 'axial mass' in GeV
static const double kQelFA0       = -1.23;   // Quasi-Elastic Axial Form Factor FA(q2=0)
static const double kQelMv2       =  TMath::Power(kQelMv,2);
static const double kQelMa2       =  TMath::Power(kQelMa,2);

//Resonance:
static const double kResMv        =  0.840;  // Resonanve 'vector mass' in GeV
static const double kResMa        =  1.032;  // Resonance 'axial mass' in GeV
static const double kResMv2       =  TMath::Power(kResMv,2);
static const double kResMa2       =  TMath::Power(kResMa,2);

//Coherent:
static const double kCohMa        =  1.032;  // Coherent 'axial mass' in GeV
static const double kCohMa2       =  TMath::Power(kCohMa,2);
static const double kCohR0        =  1.000;  // Nuclear size scale in COH model
static const double kCohReImAmpl  =  0.300;  // Re/Im of forward pion scattering amplitude

//----- constants in Feynman-Kislinger-Ravndal model for baryon excitation

static const double kOmega = 1.05;
static const double kZeta  = 0.75;

//----- constants in Rein-Seghal model for Resonant Single Pion Production

static const double kReinSeghalWCut = 1.4; // GeV (can be overriden in RS Alg)

//----- constants used in Paschos-Lalakulich v RES model for P33(1232)
//
const double kPlResMa          = 1.05; // 'axial mass' in GeV
const double kPlResMv          = 0.84; // 'vector mass' in GeV

const double kPlRes_f_pi       = 0.97*kPionMass;   // f_pi, GeV

const double kPlRes_f3_P1232_V =  1.95/kNucleonMass;
const double kPlRes_f4_P1232_V = -1.95/kNucleonMass;
const double kPlRes_f5_P1232_V =  0;
const double kPlRes_f5_P1232_A =  1.2;
const double kPlRes_f4_P1232_A = -0.3/kNucleonMass/kNucleonMass;
const double kPlRes_f3_P1232_A =  0;
const double kPlRes_f6_P1232_A =  kPlRes_f5_P1232_A;

//----- Constants used in hadronization methods

static const int kMaxMultiplicity = 15;

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

//----- Standard model constants used in Parton Model

static const double kGL       =  (1./2.) - (2./3.)* kSin8w_2;
static const double kGR       =          - (2./3.)* kSin8w_2;
static const double kGLprime  = -(1./2.) + (1./3.)* kSin8w_2;
static const double kGRprime  =            (1./3.)* kSin8w_2;
static const double kGL2      = TMath::Power(kGL,2);
static const double kGR2      = TMath::Power(kGR,2);
static const double kGLprime2 = TMath::Power(kGLprime,2);
static const double kGRprime2 = TMath::Power(kGRprime,2);

//----- Misc constants for empirical formulas

// Ro in nuclear radius formula R=Ro*A^(1/3), in GeV^-1
static const double kNucRo = 1.2E-15*m;

//----- Misc hard limits, cuts

static const double kMinQ2Limit = 1e-6;  // GeV^2

//----- Misc non-physics model constans

// maximum allowed number of iterations in rejection MC method
// before selecting a valid number
static const unsigned int kRjMaxIterations = 1000;

// Default random number generator seed number. It can be overriden
// setting the $GSEED env. var. or by using RandomGen::SetSeed(int)
static const unsigned int kDefaultRandSeed = 65539;

} // namespace constants
} // namespace genie

#endif // _CONSTANTS_H_


