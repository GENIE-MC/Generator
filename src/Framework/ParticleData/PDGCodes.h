//____________________________________________________________________________
/*!

\file     PDGCodes.h

\brief    Most commonly used PDG codes.
          A set of utility functions to handle PDG codes is provided in PDGUtils

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _PDG_CODES_H_
#define _PDG_CODES_H_

namespace genie {

const int kPdgNuE              =  12;   //
const int kPdgAntiNuE          = -12;   //
const int kPdgNuMu             =  14;   //
const int kPdgAntiNuMu         = -14;   // 
const int kPdgNuTau            =  16;   //
const int kPdgAntiNuTau        = -16;   //

const int kPdgElectron         =  11;   //
const int kPdgPositron         = -11;   //
const int kPdgMuon             =  13;   //
const int kPdgAntiMuon         = -13;   //
const int kPdgTau              =  15;   //
const int kPdgAntiTau          = -15;   //

const int kPdgUQuark           =   2;   //
const int kPdgAntiUQuark       =  -2;   //
const int kPdgDQuark           =   1;   //
const int kPdgAntiDQuark       =  -1;   //
const int kPdgSQuark           =   3;   //
const int kPdgAntiSQuark       =  -3;   //
const int kPdgCQuark           =   4;   //
const int kPdgAntiCQuark       =  -4;   //
const int kPdgBQuark           =   5;   //
const int kPdgAntiBQuark       =  -5;   //
const int kPdgTQuark           =   6;   //
const int kPdgAntiTQuark       =  -6;   //

const int kPdgDDDiquarkS1      =  1103; // dd, spin = 1 
const int kPdgUDDiquarkS0      =  2101; // ud, spin = 0 
const int kPdgUDDiquarkS1      =  2103; // ud, spin = 1 
const int kPdgUUDiquarkS1      =  2203; // uu, spin = 1 
const int kPdgSDDiquarkS0      =  3101; // sd, spin = 0
const int kPdgSDDiquarkS1      =  3103; // sd, spin = 1
const int kPdgSUDiquarkS0      =  3201; // su, spin = 0
const int kPdgSUDiquarkS1      =  3203; // su, spin = 1
const int kPdgSSDiquarkS1      =  3303; // ss, spin = 1

const int kPdgProton           =  2212; //
const int kPdgAntiProton       = -2212; //
const int kPdgNeutron          =  2112; //
const int kPdgAntiNeutron      = -2112; //
const int kPdgLambda           =  3122; // Lambda
const int kPdgAntiLambda       = -3122; // \bar{Lambda}
const int kPdgSigmaP           =  3222; // Sigma+
const int kPdgSigma0           =  3212; // Sigma0
const int kPdgSigmaM           =  3112; // Sigma-
const int kPdgAntiSigmaP       = -3222; // \bar{Sigma+}
const int kPdgAntiSigma0       = -3212; // \bar{Sigma0}
const int kPdgAntiSigmaM       = -3112; // \bar{Sigma-}
const int kPdgXi0              =  3322; // Xi0
const int kPdgXiM              =  3312; // Xi-
const int kPdgAntiXi0          = -3322; // \bar{Xi0}
const int kPdgAntiXiP          = -3312; // \bar{Xi+}
const int kPdgOmegaM           =  3334; // Omega-
const int kPdgAntiOmegaP       = -3334; // \bar{Omega+}
const int kPdgLambdaPc         =  4122; // Lambda+_{c}
const int kPdgSigma0c          =  4112; // Sigma0_{c}
const int kPdgSigmaPc          =  4212; // Sigma+_{c}
const int kPdgSigmaPPc         =  4222; // Sigma++_{c}

const int kPdgP33m1232_DeltaM  =  1114; // P33(1232) Delta-
const int kPdgP33m1232_Delta0  =  2114; // P33(1232) Delta0
const int kPdgP33m1232_DeltaP  =  2214; // P33(1232) Delta+
const int kPdgP33m1232_DeltaPP =  2224; // P33(1232) Delta++
const int kPdgS11m1535_N0      = 22112; // S11(1535) N0
const int kPdgS11m1535_NP      = 22212; // S11(1535) N+
const int kPdgD13m1520_N0      =  1214; // D13(1520) N0
const int kPdgD13m1520_NP      =  2124; // D13(1520) N+
const int kPdgS11m1650_N0      = 32112; // S11(1650) N0
const int kPdgS11m1650_NP      = 32212; // S11(1650) N+
const int kPdgD13m1700_N0      = 21214; // D13(1700) N0
const int kPdgD13m1700_NP      = 22124; // D13(1700) N+
const int kPdgD15m1675_N0      =  2116; // D15(1675) N0
const int kPdgD15m1675_NP      =  2216; // D15(1675) N+
const int kPdgS31m1620_DeltaM  =  1112; // S31(1620) Delta-
const int kPdgS31m1620_Delta0  =  1212; // S31(1620) Delta0
const int kPdgS31m1620_DeltaP  =  2122; // S31(1620) Delta+
const int kPdgS31m1620_DeltaPP =  2222; // S31(1620) Delta++
const int kPdgD33m1700_DeltaM  = 11114; // D33(1700) Delta-
const int kPdgD33m1700_Delta0  = 12114; // D33(1700) Delta0
const int kPdgD33m1700_DeltaP  = 12214; // D33(1700) Delta+
const int kPdgD33m1700_DeltaPP = 12224; // D33(1700) Delta++
const int kPdgP11m1440_N0      = 12112; // P11(1440) N0
const int kPdgP11m1440_NP      = 12212; // P11(1440) N+
const int kPdgP13m1720_N0      = 31214; // P13(1720) N0
const int kPdgP13m1720_NP      = 32124; // P13(1720) N+
const int kPdgF15m1680_N0      = 12116; // F15(1680) N0
const int kPdgF15m1680_NP      = 12216; // F15(1680) N+
const int kPdgP31m1910_DeltaM  = 21112; // P31(1910) Delta-
const int kPdgP31m1910_Delta0  = 21212; // P31(1910) Delta0
const int kPdgP31m1910_DeltaP  = 22122; // P31(1910) Delta+
const int kPdgP31m1910_DeltaPP = 22222; // P31(1910) Delta++
const int kPdgP33m1920_DeltaM  = 21114; // P33(1920) Delta-
const int kPdgP33m1920_Delta0  = 22114; // P33(1920) Delta0
const int kPdgP33m1920_DeltaP  = 22214; // P33(1920) Delta+
const int kPdgP33m1920_DeltaPP = 22224; // P33(1920) Delta++
const int kPdgF35m1905_DeltaM  =  1116; // F35(1905) Delta-
const int kPdgF35m1905_Delta0  =  1216; // F35(1905) Delta0
const int kPdgF35m1905_DeltaP  =  2126; // F35(1905) Delta+
const int kPdgF35m1905_DeltaPP =  2226; // F35(1905) Delta++
const int kPdgF37m1950_DeltaM  =  1118; // F37(1950) Delta-
const int kPdgF37m1950_Delta0  =  2118; // F37(1950) Delta0
const int kPdgF37m1950_DeltaP  =  2218; // F37(1950) Delta+
const int kPdgF37m1950_DeltaPP =  2228; // F37(1950) Delta++
const int kPdgP11m1710_N0      = 42112; // P11(1710) N0
const int kPdgP11m1710_NP      = 42212; // P11(1710) N+

const int kPdgPiP              =   211; // pi+
const int kPdgPiM              =  -211; // pi-
const int kPdgPi0              =   111; // pi0
const int kPdgEta              =   221; // eta
const int kPdgEtaPrm           =   331; // eta' (prime)
const int kPdgEtac             =   441; // eta_{c}
const int kPdgEtab             =   551; // eta_{b}
const int kPdgRhoP             =   213; // rho+
const int kPdgRhoM             =  -213; // rho-
const int kPdgRho0             =   113; // rho0
const int kPdgomega            =   223; // omega (the meson, not Omega the baryon)
const int kPdgPhi              =   333; // phi
const int kPdgJpsi             =   443; // J/psi
const int kPdgY                =   553; // Y
const int kPdgKP               =   321; // K+
const int kPdgKM               =  -321; // K-
const int kPdgK0               =   311; // K0
const int kPdgAntiK0           =  -311; // \bar{K0}
const int kPdgK0L              =   130; // K0_{long}
const int kPdgK0S              =   310; // K0_{short}
const int kPdgKStarP           =   323; // Kstar+(892)
const int kPdgKStarM           =  -323; // Kstar-(892) 
const int kPdgKStar0           =   313; // Kstar0(892)
const int kPdgDP               =   411; // D+
const int kPdgDM               =  -411; // D-
const int kPdgD0               =   421; // D0
const int kPdgAntiD0           =  -421; // \bar{D0}
const int kPdgDPs              =   431; // D+_{s}
const int kPdgDMs              =  -431; // D-_{s}

const int kPdgGluon            =    21; // gluon
const int kPdgGamma            =    22; // photon
const int kPdgZ0               =    23; // Z
const int kPdgWP               =    24; // W+
const int kPdgWM               =   -24; // W-

//
// Note: PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z)
// PDG2006 convention: 10LZZZAAAI 
// Define names for some commonly used nuclear PDG codes:
//
const int kPdgTgtFreeP     = 1000010010;
const int kPdgTgtFreeN     = 1000000010;
const int kPdgTgtDeuterium = 1000010020;
const int kPdgTgtC12       = 1000060120;
const int kPdgTgtO16       = 1000080160;
const int kPdgTgtCa40      = 1000200400;
const int kPdgTgtFe56      = 1000260560;

//
// PDG codes for GENIE special particles
//
const int kPdgHadronicSyst    = 2000000001; // DIS hadronic system before hadronization
const int kPdgHadronicBlob    = 2000000002; // Unmodelled fraction of the hadronic system
const int kPdgBindino         = 2000000101; // Binding energy subtracted from f/s nucleons
const int kPdgCoulobtron      = 2000000102; // Coulomb energy subtracted from f/s leptons
const int kPdgClusterNN       = 2000000200; // A nn cluster within a nucleus
const int kPdgClusterNP       = 2000000201; // A np cluster within a nucleus
const int kPdgClusterPP       = 2000000202; // A pp cluster within a nucleus
const int kPdgCompNuclCluster = 2000000300; // Nucleon cluster before phase decay 
const int kPdgDarkMatter      = 2000010000; // Dark matter particle for GENIE Boosted Dark Matter mode
const int kPdgMediator        = 2000010001; // Mediator particle for GENIE Boosted Dark Matter mode

//
// PDG codes for special particles used by external generators interfaced with GENIE
//
const int kPdgCluster      = 91; // PYTHIA cluster pseudo-particle
const int kPdgString       = 92; // PYTHIA string pseudo-particle 
const int kPdgIndep        = 93; // PYTHIA independent fragmentation pseudo-particle 

}      // genie namespace

#endif // _PDG_CODES_H_
