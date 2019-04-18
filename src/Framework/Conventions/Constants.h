//____________________________________________________________________________
/*!

\namespace genie::constants

\brief     Basic constants

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

\created   May 03, 2004

\cpright   Copyright (c) 2003-2019, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <TMath.h>

#include "Framework/Conventions/Units.h"

namespace genie {
namespace constants {

//
// Fundamental constants
//
static const double kLightSpeed    = 1.;
static const double kPlankConstant = 1.;

//
// pi, e,...
//
static const double kPi    = TMath::Pi();
static const double kPi2   = TMath::Power(kPi,2);
static const double kPi3   = TMath::Power(kPi,3);
static const double kPi4   = TMath::Power(kPi,4);
static const double kSqrtPi= TMath::Sqrt(kPi);

static const double kNapierConst     = TMath::E();
static const double kSqrtNapierConst = TMath::Sqrt(kNapierConst);

//
// Avogadro number, compton wavelength and such...
//
static const double kNA    = 6.02214179E+23;             
static const double kLe    = 3.8616E-11 *units::cm;
static const double kLe2   = TMath::Power(kLe,2);          

//
// Coupling constants
//
static const double kAem   = 1./137.03599976;       // EM coupling const, dimensionless
static const double kAem2  = TMath::Power(kAem,2);
static const double kGF    = 1.16639E-5;            // Fermi const from b-decay, in GeV^-2
static const double kGF2   = TMath::Power(kGF,2);

//
// Masses
//
// For simplicity, the most commonly used particle masses defined here.
// In general, however, particle masses in GENIE classes should be obtained
// through the genie::PDGLibrary as shown below:
// double mass = PDGLibrary::Instance()->Find(pdg_code)->Mass();
// For consistency, the values below must match whatever is used in PDGLibrary.
//
static const double kElectronMass   =  5.109989461e-04;        // GeV
static const double kMuonMass       =  1.056583745e-01;         // GeV
static const double kTauMass        =  1.77686e+00;             // GeV
static const double kPionMass       =  1.3957018e-01;          // GeV
static const double kPi0Mass        =  1.349766e-01;           // GeV
static const double kProtonMass     =  9.38272081e-01;           // GeV
static const double kNeutronMass    =  9.39565413e-01;           // GeV
static const double kNucleonMass    =  (kProtonMass+kNeutronMass)/2.;  //GeV
static const double kLightestChmHad =  1.870;               // GeV ~lightest charm hadron+
static const double kPhotontest     =  1E-6;                // GeV



static const double kElectronMass2  =  TMath::Power(kElectronMass,2); // GeV^2
static const double kMuonMass2      =  TMath::Power(kMuonMass,2);     // GeV^2
static const double kTauMass2       =  TMath::Power(kTauMass,2);      // GeV^2
static const double kPionMass2      =  TMath::Power(kPionMass,2);     // GeV^2
static const double kProtonMass2    =  TMath::Power(kProtonMass,2);   // GeV^2
static const double kNeutronMass2   =  TMath::Power(kNeutronMass,2);  // GeV^2
static const double kNucleonMass2   =  TMath::Power(kNucleonMass,2);  // GeV^2

static const double kMw             =  8.0385e+01;                // GeV - W boson mass
static const double kMz             =  9.11876e+01;                // GeV - Z boson mass
static const double kMw2            =  TMath::Power(kMw,2);  // GeV^2
static const double kMz2            =  TMath::Power(kMz,2);  // GeV^2

//
// Misc constants for empirical formulas
//
static const double kNucRo      = 1.2E-15 * units::m;            // Ro in nuclear radius formula R=Ro*A^(1/3), in GeV^-1
static const double kNucDensity = 2.3E+17 * units::kg/units::m3; // Nuclear density (in nuclear core), in GeV^4


//FMTOGEV
static const double FMTOGEV= 5.0761421;


//
// Earth consts
//
static const double kREarth = 6371 * units::km; 

//
// Sqrts frequently encountered in helicity amplitude calculations
//
static const double kSqrt2     =  1.4142136;
static const double kSqrt3     =  1.7320508;
static const double kSqrt4     =  2.0;
static const double kSqrt5     =  2.236068;
static const double kSqrt6     =  2.4494897;
static const double kSqrt7     =  2.6457513;
static const double kSqrt8     =  2.8284271;
static const double kSqrt9     =  3.0;
static const double kSqrt10    =  3.1622777;
static const double kSqrt12    =  3.4641016;
static const double kSqrt15    =  3.8729833;
static const double kSqrt18    =  4.2426407;
static const double kSqrt20    =  4.472136;
static const double kSqrt24    =  4.8989795;
static const double kSqrt27    =  5.1961524;
static const double kSqrt30    =  5.4772256;
static const double kSqrt35    =  5.9160798;
static const double kSqrt40    =  6.3245553;
static const double kSqrt60    =  7.7459667;
static const double kSqrt120   = 10.954451;
static const double k1_Sqrt2   =  0.70710678;
static const double k1_Sqrt3   =  0.57735027;
static const double k1_Sqrt5   =  0.44721360;
static const double k1_Sqrt6   =  0.40824829;
static const double k1_Sqrt7   =  0.37796447;
static const double k1_Sqrt10  =  0.31622777;
static const double k1_Sqrt15  =  0.25819889;
static const double k1_Sqrt24  =  0.20412415;
static const double k1_Sqrt30  =  0.18257419;
static const double k1_Sqrt35  =  0.16903085;
static const double k1_Sqrt60  =  0.12909944;
static const double k1_Sqrt120 =  0.091287093;
static const double k2_Sqrt3   =  1.1547005;
static const double k2_Sqrt5   =  0.89442719;
static const double k2_Sqrt15  =  0.51639778;
static const double k2_Sqrt35  =  0.3380617;
static const double k3_Sqrt2   =  2.1213203;
static const double k3_Sqrt5   =  1.3416408;
static const double k3_Sqrt10  =  0.9486833;
static const double k3_Sqrt20  =  0.67082039;
static const double k3_Sqrt40  =  0.47434165;
static const double kSqrt2_3   =  0.81649658;
static const double kSqrt2_5   =  0.63245553;
static const double kSqrt2_6   =  0.57735027;
static const double kSqrt2_7   =  0.53452248;
static const double kSqrt2_15  =  0.36514837;
static const double kSqrt3_2   =  1.2247449;
static const double kSqrt3_4   =  0.8660254;
static const double kSqrt3_5   =  0.77459667;
static const double kSqrt3_8   =  0.61237244;
static const double kSqrt3_10  =  0.54772256;
static const double kSqrt3_18  =  0.40824829;
static const double kSqrt3_20  =  0.38729833;
static const double kSqrt3_35  =  0.29277002;
static const double kSqrt3_40  =  0.27386128;
static const double kSqrt4_15  =  0.51639778;
static const double kSqrt5_2   =  1.5811388;
static const double kSqrt5_3   =  1.2909944;
static const double kSqrt5_8   =  0.79056942;
static const double kSqrt5_12  =  0.64549722;
static const double kSqrt6_5   =  1.0954451;
static const double kSqrt6_35  =  0.41403934;
static const double kSqrt9_10  =  0.9486833;
static const double kSqrt9_40  =  0.47434165;
static const double kSqrt18_5  =  1.8973666;
static const double kSqrt18_20 =  0.9486833;
static const double kSqrt18_35 =  0.71713717;
static const double kSqrt24_35 =  0.82807867;
static const double kSqrt27_10 =  1.6431677;
static const double kSqrt27_40 =  0.82158384;

} // namespace constants
} // namespace genie

#endif // _CONSTANTS_H_


