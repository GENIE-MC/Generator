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

/*
PDG naming convention:
-----------------------
Mesons are represented by the form NML and baryons by NMKL where N, M and K 
are the quark numbers (1=d, 2=u, 3=s, 4=c, 5=b, 6=t) and L=2J+1 where J is 
the spin. K short and K long are exceptions being 310 and 130. The quarks 
are in decreasing mass order from left to right except that Lambda (3122) 
has the order switched to distinguish it from Sigma (3212). The fifth digit 
is used to distinguish particles with the same spin and quark content. 
Particles have positive numbers, antiparticles have the negative of their 
counterparts. The positive kaon and positive B meson are particles. 
The electron, electron neutrino, negative muon etc are 11, 12, 13 etc. 
The gluon, photon, Z, and W are 21, 22, 23 and 24. 
The PDG naming convention does not cover nuclei. GENIE uses the MINOS
naming convention to code nuclei.

MINOS naming conventions (for nuclei):
---------------------------------------
1AAAZZZ000
*/

const int kPdgNuE          =  12;
const int kPdgAntiNuE      = -12;
const int kPdgNuMu         =  14;
const int kPdgAntiNuMu     = -14;
const int kPdgNuTau        =  16;
const int kPdgAntiNuTau    = -16;

const int kPdgElectron     =  11;
const int kPdgPositron     = -11;
const int kPdgMuon         =  13;
const int kPdgAntiMuon     = -13;
const int kPdgTau          =  15;
const int kPdgAntiTau      = -15;

const int kPdgUQuark       =   2;
const int kPdgAntiUQuark   =  -2;
const int kPdgDQuark       =   1;  
const int kPdgAntiDQuark   =  -1;
const int kPdgSQuark       =   3;
const int kPdgAntiSQuark   =  -3;
const int kPdgCQuark       =   4;
const int kPdgAntiCQuark   =  -4;
const int kPdgBQuark       =   5;
const int kPdgAntiBQuark   =  -5;
const int kPdgTQuark       =   6;
const int kPdgAntiTQuark   =  -6;
const int kPdgUUDiquarkS1  =  2203; // uu, spin = 1 - triplet
const int kPdgUDDiquarkS0  =  2101; // ud, spin = 0 - singlet
const int kPdgUDDiquarkS1  =  2103; // ud, spin = 1 - triplet
const int kPdgDDDiquarkS1  =  1103; // dd, spin = 1 - triplet

const int kPdgProton       =  2212;
const int kPdgAntiProton   = -2212;
const int kPdgNeutron      =  2112;
const int kPdgAntiNeutron  = -2112;
const int kPdgLambda       =  3122; // Lambda
const int kPdgAntiLambda   = -3122; // \bar{Lambda}
const int kPdgSigmaP       =  3222; // Sigma+
const int kPdgSigma0       =  3212; // Sigma0
const int kPdgSigmaM       =  3112; // Sigma-
const int kPdgAntiSigmaP   = -3112; // Sigma+
const int kPdgAntiSigma0   = -3212; // Sigma0
const int kPdgAntiSigmaM   = -3112; // Sigma-
const int kPdgXi0          =  3322; // Xi0
const int kPdgXiM          =  3312; // Xi-
const int kPdgAntiXi0      = -3322; // \bar{Xi0}
const int kPdgAntiXiP      = -3312; // \bar{Xi+}
const int kPdgOmegaM       =  3332; // Omega-
const int kPdgAntiOmegaP   = -3332; // \bar{Omega+}
const int kPdgLambdaPc     =  4122; // Lambda+_{c}
const int kPdgSigmaPc      =  4212; // Sigma+_{c}
const int kPdgSigmaPPc     =  4222; // Sigma++_{c}

const int kPdgPiP          =   211; // pi+
const int kPdgPiM          =  -211; // pi-
const int kPdgPi0          =   111; // pi0
const int kPdgEta          =   221; // eta
const int kPdgEtaPrm       =   331; // eta' (prime)
const int kPdgEtac         =   441; // eta_{c}
const int kPdgEtab         =   551; // eta_{b}
const int kPdgRhoP         =   213; // rho+
const int kPdgRhoM         =  -213; // rho-
const int kPdgRho0         =   113; // rho0
const int kPdgomega        =   223; // omega (the meson, not Omega the baryon)
const int kPdgPhi          =   333; // phi
const int kPdgJpsi         =   443; // J/psi
const int kPdgY            =   553; // Y
const int kPdgKP           =   321; // K+
const int kPdgKM           =  -321; // K-
const int kPdgK0           =   311; // K0
const int kPdgK0L          =   130; // K0_{long}
const int kPdgK0S          =   310; // K0_{short}
const int kPdgDP           =   411; // D+
const int kPdgDM           =  -411; // D-
const int kPdgD0           =   421; // D0
const int kPdgAntiD0       =   421; // \bar{D0} - no separate PDG code?
const int kPdgDPs          =   431; // D+_{s}
const int kPdgDMs          =  -431; // D-_{s}

const int kPdgGluon        =    21; // gluon
const int kPdgGamma        =    22; // photon
const int kPdgZ0           =    23; // Z
const int kPdgWP           =    24; // W+
const int kPdgWM           =   -24; // W-

// note: PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z);
const int kPdgTgtFreeP     = 1001001000;
const int kPdgTgtFreeN     = 1001000000;
const int kPdgTgtFe56      = 1056026000;
const int kPdgTgtC12       = 1012006000;
const int kPdgTgtO16       = 1016008000;

// GENIE special particles
const int kPdgBindino      = 1111111001; 
const int kPdgHadronicSyst = 1111111002; 
const int kPdgHadronicBlob = 1111111003; // not explicitly hadronized fraction of hadronic system

// PYTHIA special particles
const int kPdgCluster      = 91; 
const int kPdgString       = 92; 
const int kPdgIndep        = 93; 


}      // genie namespace

#endif // _PDG_CODES_H_
