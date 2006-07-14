//____________________________________________________________________________
/*!

\file     PDGCodes.h

\brief    Most commonly used PDG codes.
          A set of utility functions to hagle PDG codes is provided in
          PDGUtils

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _PDG_CODES_H_
#define _PDG_CODES_H_

namespace genie {

const int kPdgNuE          =  12;
const int kPdgNuEBar       = -12;
const int kPdgNuMu         =  14;
const int kPdgNuMuBar      = -14;
const int kPdgNuTau        =  16;
const int kPdgNuTauBar     = -16;
const int kPdgElectron     =  11;
const int kPdgPositron     = -11;
const int kPdgMuon         =  13;
const int kPdgAntiMuon     = -13;
const int kPdgTau          =  15;
const int kPdgAntiTau      = -15;
const int kPdgUQuark       =   2;
const int kPdgUQuarkBar    =  -2;
const int kPdgDQuark       =   1;  
const int kPdgDQuarkBar    =  -1;
const int kPdgSQuark       =   3;
const int kPdgSQuarkBar    =  -3;
const int kPdgCQuark       =   4;
const int kPdgCQuarkBar    =  -4;
const int kPdgBQuark       =   5;
const int kPdgBQuarkBar    =  -5;
const int kPdgTQuark       =   6;
const int kPdgTQuarkBar    =  -6;
const int kPdgUUDiquarkS1  =  2203; // uu, spin = 1 - triplet
const int kPdgUDDiquarkS0  =  2101; // ud, spin = 0 - singlet
const int kPdgUDDiquarkS1  =  2103; // ud, spin = 1 - triplet
const int kPdgDDDiquarkS1  =  1103; // dd, spin = 1 - triplet
const int kPdgProton       =  2212;
const int kPdgNeutron      =  2112;
const int kPdgPiPlus       =   211;
const int kPdgPiMinus      =  -211;
const int kPdgPi0          =   111;
const int kPdgKPlus        =   321;
const int kPdgKMinus       =  -321;
const int kPdgK0           =   311;
const int kPdgLambdacP     =  4122; // #Lambda_{c}^{+}
const int kPdgSigmacP      =  4212; // #Sigma_{c}^{+}
const int kPdgSigmacPP     =  4222; // #Sigma_{c}^{++}
const int kPdgDPlus        =   411; // D^{+}
const int kPdgDMinus       =  -411; // D^{-}
const int kPdgD0           =   421; // D^{0}
const int kPdgD0Bar        =   421; // #bar{D^{0}} - no separate PDG code?
const int kPdgDsPlus       =   431; // D_{s}^{+}
const int kPdgDsMinus      =  -431; // D_{s}^{-}
const int kPdgGamma        =    22;
const int kPdgZ0           =    23;
const int kPdgWPlus        =    24;
const int kPdgWMinus       =   -24;

// note: PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z);
const int kPdgTgtFreeP     = 1001001000;
const int kPdgTgtFreeN     = 1001000000;
const int kPdgTgtFe56      = 1056026000;
const int kPdgTgtC12       = 1012006000;
const int kPdgTgtO16       = 1016008000;

// GENIE special particles
const int kPdgBindino      = 1111111001; 
const int kPdgHadronicSyst = 1111111002; 

// PYTHIA special particles
const int kPdgCluster      = 91; 
const int kPdgString       = 92; 
const int kPdgIndep        = 93; 


}      // genie namespace

#endif // _PDG_CODES_H_
