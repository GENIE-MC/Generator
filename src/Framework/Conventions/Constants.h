//____________________________________________________________________________
/*!

\namespace genie::constants

\brief     Basic constants

\author    Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   May 03, 2004

\cpright   Copyright (c) 2003-2020, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org     
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
static const double kAem   = 1./137.035999139;       // EM coupling const, dimensionless
static const double kAem2  = TMath::Power(kAem,2);
static const double kGF    =  1.1663787E-5;            // Fermi const from b-decay, in GeV^-2
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
static const double kMuonMass       =  1.056583745e-01;        // GeV
static const double kTauMass        =  1.77686e+00;            // GeV
static const double kPionMass       =  1.3957018e-01;          // GeV
static const double kPi0Mass        =  1.349766e-01;           // GeV
static const double kProtonMass     =  9.38272081e-01;         // GeV
static const double kNeutronMass    =  9.39565413e-01;         // GeV
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

static const double kMw             =  8.0379e+01;           // GeV - W boson mass
static const double kMz             =  9.11876e+01;          // GeV - Z boson mass
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
static const double kSqrt2     =  TMath::Sqrt2();
static const double kSqrt3     =  TMath::Sqrt(3);
static const double kSqrt4     =  TMath::Sqrt(4);
static const double kSqrt5     =  TMath::Sqrt(5);
static const double kSqrt6     =  TMath::Sqrt(6);
static const double kSqrt7     =  TMath::Sqrt(7);
static const double kSqrt8     =  TMath::Sqrt(8);
static const double kSqrt9     =  TMath::Sqrt(9);
static const double kSqrt10    =  TMath::Sqrt(10);
static const double kSqrt12    =  TMath::Sqrt(12);
static const double kSqrt15    =  TMath::Sqrt(15);
static const double kSqrt18    =  TMath::Sqrt(18);
static const double kSqrt20    =  TMath::Sqrt(20);
static const double kSqrt24    =  TMath::Sqrt(24);
static const double kSqrt27    =  TMath::Sqrt(27);
static const double kSqrt30    =  TMath::Sqrt(30);
static const double kSqrt35    =  TMath::Sqrt(35);
static const double kSqrt40    =  TMath::Sqrt(40);
static const double kSqrt60    =  TMath::Sqrt(60);
static const double kSqrt120   =  TMath::Sqrt(120);
static const double k1_Sqrt2   =  1/TMath::Sqrt2();
static const double k1_Sqrt3   =  1/TMath::Sqrt(3);
static const double k1_Sqrt5   =  1/TMath::Sqrt(5);
static const double k1_Sqrt6   =  1/TMath::Sqrt(6);
static const double k1_Sqrt7   =  1/TMath::Sqrt(7);
static const double k1_Sqrt10  =  1/TMath::Sqrt(10);
static const double k1_Sqrt15  =  1/TMath::Sqrt(15);
static const double k1_Sqrt24  =  1/TMath::Sqrt(24);
static const double k1_Sqrt30  =  1/TMath::Sqrt(30);
static const double k1_Sqrt35  =  1/TMath::Sqrt(35);
static const double k1_Sqrt60  =  1/TMath::Sqrt(60);
static const double k1_Sqrt120 =  1/TMath::Sqrt(120);
static const double k2_Sqrt3   =  2/TMath::Sqrt(3);
static const double k2_Sqrt5   =  2/TMath::Sqrt(5);
static const double k2_Sqrt15  =  2/TMath::Sqrt(15);
static const double k2_Sqrt35  =  2/TMath::Sqrt(35);
static const double k3_Sqrt2   =  3/TMath::Sqrt2();
static const double k3_Sqrt5   =  3/TMath::Sqrt(5);
static const double k3_Sqrt10  =  3/TMath::Sqrt(10);
static const double k3_Sqrt20  =  3/TMath::Sqrt(20);
static const double k3_Sqrt40  =  3/TMath::Sqrt(40);
static const double kSqrt2_3   =  TMath::Sqrt(2.0/3);
static const double kSqrt2_5   =  TMath::Sqrt(2.0/5);
static const double kSqrt2_6   =  TMath::Sqrt(2.0/6);
static const double kSqrt2_7   =  TMath::Sqrt(2.0/7);
static const double kSqrt2_15  =  TMath::Sqrt(2.0/15);
static const double kSqrt3_2   =  TMath::Sqrt(3.0/2);
static const double kSqrt3_4   =  TMath::Sqrt(3.0/4);
static const double kSqrt3_5   =  TMath::Sqrt(3.0/5);
static const double kSqrt3_8   =  TMath::Sqrt(3.0/8);
static const double kSqrt3_10  =  TMath::Sqrt(3.0/10);
static const double kSqrt3_18  =  TMath::Sqrt(3.0/18);
static const double kSqrt3_20  =  TMath::Sqrt(3.0/20);
static const double kSqrt3_35  =  TMath::Sqrt(3.0/35);
static const double kSqrt3_40  =  TMath::Sqrt(3.0/40);
static const double kSqrt4_15  =  TMath::Sqrt(4.0/15);
static const double kSqrt5_2   =  TMath::Sqrt(5.0/2);
static const double kSqrt5_3   =  TMath::Sqrt(5.0/3);
static const double kSqrt5_8   =  TMath::Sqrt(5.0/8);
static const double kSqrt5_12  =  TMath::Sqrt(5.0/12);
static const double kSqrt6_5   =  TMath::Sqrt(6.0/5);
static const double kSqrt6_35  =  TMath::Sqrt(6.0/35);
static const double kSqrt9_10  =  TMath::Sqrt(9.0/10);
static const double kSqrt9_40  =  TMath::Sqrt(9.0/40);
static const double kSqrt18_5  =  TMath::Sqrt(18.0/5);
static const double kSqrt18_20 =  TMath::Sqrt(18.0/20);
static const double kSqrt18_35 =  TMath::Sqrt(18.0/35);
static const double kSqrt24_35 =  TMath::Sqrt(24.0/35);
static const double kSqrt27_10 =  TMath::Sqrt(27.0/10);
static const double kSqrt27_40 =  TMath::Sqrt(27.0/40);

} // namespace constants
} // namespace genie

#endif // _CONSTANTS_H_
