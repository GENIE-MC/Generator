//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <cassert>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

//____________________________________________________________________________
bool genie::pdg::IsPseudoParticle(int pdgc)
{
// ROOT's rootino has PDG code=0
// GENIE pseudoparticles are in the 2000000000-2000100000 range
// Include PYTHIA's pseudoparticles

  bool is_fake =
     ( (pdgc == 0) ||
       (pdgc > 2000000000 && pdgc < 2000100000) ||
       (pdgc == kPdgCluster || pdgc == kPdgString || pdgc == kPdgIndep)
     );
      
  return is_fake;
}
//____________________________________________________________________________
bool genie::pdg::IsIon(int pdgc)
{
  return (pdgc > 1000000000 && pdgc < 1999999999);
}
//____________________________________________________________________________
bool genie::pdg::IsParticle(int pdgc)
{
  if( genie::pdg::IsPseudoParticle (pdgc) ) return false;
  if( genie::pdg::IsIon            (pdgc) ) return false;

  return true;
}
//____________________________________________________________________________
int genie::pdg::IonPdgCodeToZ(int ion_pdgc)
{
// Decoding Z from the PDG code (PDG ion code convention: 10LZZZAAAI)

  int Z = (ion_pdgc/10000) - 1000*(ion_pdgc/10000000); // don't factor out!
  return Z;
}
//____________________________________________________________________________
int genie::pdg::IonPdgCodeToA(int ion_pdgc)
{
// Decoding A from the PDG code (PDG ion code convention: 10LZZZAAAI)

  int A = (ion_pdgc/10) - 1000*(ion_pdgc/10000); // don't factor out!
  return A;
}
//____________________________________________________________________________
int genie::pdg::IonPdgCode(int A, int Z)
{
// Build ion PDG code from A,Z

  return IonPdgCode(A,Z,0,0);
}
//____________________________________________________________________________
int genie::pdg::IonPdgCode(int A, int Z, int L, int I)
{
// Build ion PDG code from A,Z,L,I

  int ion_pdgc = 1000000000 +  L*100000000 + Z*10000 + A*10 + I;
  return ion_pdgc;
}
//____________________________________________________________________________
bool genie::pdg::IsLepton(int pdgc)
{
  bool is_neutral_lepton = genie::pdg::IsNeutralLepton(pdgc);
  bool is_charged_lepton = genie::pdg::IsChargedLepton(pdgc);

  bool is_lepton = (is_neutral_lepton || is_charged_lepton);
  return is_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsNeutralLepton(int pdgc)
{
  bool is_neutral_lepton = IsNeutrino(pdgc) || IsAntiNeutrino(pdgc);
  return is_neutral_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsChargedLepton(int pdgc)
{
  bool is_neg_lepton = genie::pdg::IsNegChargedLepton(pdgc);
  bool is_pos_lepton = genie::pdg::IsPosChargedLepton(pdgc);

  bool is_charged_lepton = is_neg_lepton || is_pos_lepton;
  return is_charged_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsNeutrino(int pdgc)
{
  bool is_nu = (pdgc == kPdgNuE)  ||
               (pdgc == kPdgNuMu) ||
               (pdgc == kPdgNuTau);
  return is_nu;
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNeutrino(int pdgc)
{
  bool is_nubar = (pdgc == kPdgAntiNuE)  ||
                  (pdgc == kPdgAntiNuMu) ||
                  (pdgc == kPdgAntiNuTau);

  return is_nubar;
}
//____________________________________________________________________________
bool genie::pdg::IsDarkMatter(int pdgc)
{
  bool is_dm = (pdgc == kPdgDarkMatter);
  return is_dm;
}
//____________________________________________________________________________
bool genie::pdg::IsNegChargedLepton(int pdgc)
{
  bool is_neg_lepton = (pdgc ==  kPdgElectron) ||
                       (pdgc ==  kPdgMuon)     ||
                       (pdgc ==  kPdgTau);

  return is_neg_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsPosChargedLepton(int pdgc)
{
  bool is_pos_lepton = (pdgc ==  kPdgPositron) ||
                       (pdgc ==  kPdgAntiMuon) ||
                       (pdgc ==  kPdgAntiTau);

  return is_pos_lepton;
}
//____________________________________________________________________________

bool genie::pdg::IsNuE(int pdgc)
{
  return (pdgc == kPdgNuE);
}
//____________________________________________________________________________
bool genie::pdg::IsNuMu(int pdgc)
{
  return (pdgc == kPdgNuMu);
}
//____________________________________________________________________________
bool genie::pdg::IsNuTau(int pdgc)
{
  return (pdgc == kPdgNuTau);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuE(int pdgc)
{
  return (pdgc == kPdgAntiNuE);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuMu(int pdgc)
{
  return (pdgc == kPdgAntiNuMu);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuTau(int pdgc)
{
  return (pdgc == kPdgAntiNuTau);
}
//____________________________________________________________________________
bool genie::pdg::IsElectron(int pdgc)
{
  return (pdgc == kPdgElectron);
}
//____________________________________________________________________________
bool genie::pdg::IsPositron(int pdgc)
{
  return (pdgc == kPdgPositron);
}
//____________________________________________________________________________
bool genie::pdg::IsMuon(int pdgc)
{
  return (pdgc == kPdgMuon);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiMuon(int pdgc)
{
  return (pdgc == kPdgAntiMuon);
}
//____________________________________________________________________________
bool genie::pdg::IsTau(int pdgc)
{
  return (pdgc == kPdgTau);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiTau(int pdgc)
{
  return (pdgc == kPdgAntiTau);
}
//____________________________________________________________________________
int genie::pdg::Neutrino2ChargedLepton(int pdgc)
{
  switch(pdgc) {
       case (kPdgNuE)      : return kPdgElectron; break;
       case (kPdgAntiNuE)  : return kPdgPositron; break;
       case (kPdgNuMu)     : return kPdgMuon;     break;
       case (kPdgAntiNuMu) : return kPdgAntiMuon; break;
       case (kPdgNuTau)    : return kPdgTau;      break;
       case (kPdgAntiNuTau): return kPdgAntiTau;  break;
  }
  return -1;
}
//____________________________________________________________________________
bool genie::pdg::IsDiQuark(int pdgc)
{
  return ( pdgc == kPdgDDDiquarkS1 || pdgc == kPdgUDDiquarkS0 ||
           pdgc == kPdgUDDiquarkS1 || pdgc == kPdgUUDiquarkS1 ||
           pdgc == kPdgSDDiquarkS0 || pdgc == kPdgSDDiquarkS1 ||
           pdgc == kPdgSUDiquarkS0 || pdgc == kPdgSUDiquarkS1 ||
           pdgc == kPdgSSDiquarkS1
         );
}
//____________________________________________________________________________
bool genie::pdg::IsQuark(int pdgc)
{
  return ( pdgc == kPdgDQuark || pdgc == kPdgUQuark ||
           pdgc == kPdgSQuark || pdgc == kPdgCQuark ||
           pdgc == kPdgBQuark || pdgc == kPdgTQuark
         );
}
//____________________________________________________________________________
bool genie::pdg::IsAntiQuark(int pdgc)
{
  return ( pdgc == kPdgAntiDQuark || pdgc == kPdgAntiUQuark ||
           pdgc == kPdgAntiSQuark || pdgc == kPdgAntiCQuark ||
           pdgc == kPdgAntiBQuark || pdgc == kPdgAntiTQuark
         );
}
//____________________________________________________________________________
bool genie::pdg::IsUQuark(int pdgc)
{
  return (pdgc == kPdgUQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsDQuark(int pdgc)
{
  return (pdgc == kPdgDQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsSQuark(int pdgc)
{
  return (pdgc == kPdgSQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsCQuark(int pdgc)
{
  return (pdgc == kPdgCQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiUQuark(int pdgc)
{
  return (pdgc == kPdgAntiUQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiDQuark(int pdgc)
{
  return (pdgc == kPdgAntiDQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiSQuark(int pdgc)
{
  return (pdgc == kPdgAntiSQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiCQuark(int pdgc)
{
  return (pdgc == kPdgAntiCQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsPion(int pdgc)
{
  return (pdgc == kPdgPiP || pdgc == kPdgPi0 || pdgc == kPdgPiM);
}
//____________________________________________________________________________
bool genie::pdg::IsKaon(int pdgc)
{
  return (pdgc == kPdgKP || pdgc == kPdgK0 || pdgc == kPdgKM);
}
//____________________________________________________________________________
bool genie::pdg::IsProton(int pdgc)
{
  return (pdgc == kPdgProton);
}
//____________________________________________________________________________
bool genie::pdg::IsNeutron(int pdgc)
{
  return (pdgc == kPdgNeutron);
}
//____________________________________________________________________________
bool genie::pdg::IsNucleon(int pdgc)
{
  return (pdgc == kPdgNeutron || pdgc == kPdgProton);
}
//____________________________________________________________________________
bool genie::pdg::IsNeutronOrProton(int pdgc)
{
  return (pdgc == kPdgNeutron || pdgc == kPdgProton);
}
//____________________________________________________________________________
int genie::pdg::SwitchProtonNeutron(int pdgc)
{
  assert(IsProton(pdgc) || IsNeutron(pdgc));

  if (IsProton(pdgc)) return kPdgNeutron;
  else                return kPdgProton;
}
//____________________________________________________________________________
int genie::pdg::ModifyNucleonCluster(int pdgc, int dQ)
{
  assert(pdg::Is2NucleonCluster(pdgc));

  if(pdgc == kPdgClusterNN) {
    if      (dQ ==  0) { return kPdgClusterNN; }
    else if (dQ == +1) { return kPdgClusterNP; }
    else if (dQ == +2) { return kPdgClusterPP; }
    else               { return 0;             }
  }
  else
  if(pdgc == kPdgClusterNP) {
    if      (dQ == -1) { return kPdgClusterNN; }
    else if (dQ ==  0) { return kPdgClusterNP; }
    else if (dQ == +1) { return kPdgClusterPP; }
    else               { return 0;             }
  }
  else
  if(pdgc == kPdgClusterPP) {
    if      (dQ == -2) { return kPdgClusterNN; }
    else if (dQ == -1) { return kPdgClusterNP; }
    else if (dQ ==  0) { return kPdgClusterPP; }
    else               { return 0;             }
  }

  return 0;
}
//____________________________________________________________________________
bool genie::pdg::IsHadron(int pdgc)
{
  return ((pdgc>=100 && pdgc<=9999) || (pdgc>=-9999 && pdgc<=-100));
}
//____________________________________________________________________________
bool genie::pdg::IsBaryonResonance(int pdgc)
{
  return utils::res::IsBaryonResonance(pdgc);
}
//____________________________________________________________________________
bool genie::pdg::Is2NucleonCluster(int pdgc)
{
   return (
      pdgc == kPdgClusterNN   ||
      pdgc == kPdgClusterNP   ||
      pdgc == kPdgClusterPP
   );
}
//____________________________________________________________________________
int genie::pdg::GeantToPdg(int geant_code)
{
  if(geant_code ==  3) return kPdgElectron;     //    11 / e-
  if(geant_code ==  2) return kPdgPositron;     //   -11 / e+
  if(geant_code ==  6) return kPdgMuon;         //    13 / mu-
  if(geant_code ==  5) return kPdgAntiMuon;     //   -13 / mu+             
  if(geant_code == 34) return kPdgTau;          //    15 / tau-
  if(geant_code == 33) return kPdgAntiTau;      //   -15 / tau+              
  if(geant_code ==  8) return kPdgPiP;          //   211 / pi+
  if(geant_code ==  9) return kPdgPiM;          //  -211 / pi-
  if(geant_code ==  7) return kPdgPi0;          //   111 / pi0
  if(geant_code == 17) return kPdgEta;          //   221 / eta
  if(geant_code == 11) return kPdgKP;           //   321 / K+
  if(geant_code == 12) return kPdgKM;           //  -321 / K-
  if(geant_code == 10) return kPdgK0L;          //   130 / K0_{long}
  if(geant_code == 16) return kPdgK0S;          //   310 / K0_{short}
  if(geant_code == 35) return kPdgDP;           //   411 / D+
  if(geant_code == 36) return kPdgDM;           //  -411 / D-
  if(geant_code == 37) return kPdgD0;           //   421 / D0
  if(geant_code == 38) return kPdgAntiD0;       //  -421 / \bar{D0}
  if(geant_code == 39) return kPdgDPs;          //   431 / D+_{s}
  if(geant_code == 40) return kPdgDMs;          //  -431 / D-_{s}
  if(geant_code ==  1) return kPdgGamma;        //    22 / photon
  if(geant_code == 44) return kPdgZ0;           //    23 / Z
  if(geant_code == 42) return kPdgWP;           //    24 / W+
  if(geant_code == 43) return kPdgWM;           //   -24 / W-
  if(geant_code == 14) return kPdgProton;       //  2212
  if(geant_code == 15) return kPdgAntiProton;   // -2212
  if(geant_code == 13) return kPdgNeutron;      //  2112
  if(geant_code == 25) return kPdgAntiNeutron;  // -2112
  if(geant_code == 18) return kPdgLambda;       //  3122 / Lambda
  if(geant_code == 26) return kPdgAntiLambda;   // -3122 / \bar{Lambda}
  if(geant_code == 19) return kPdgSigmaP;       //  3222 / Sigma+
  if(geant_code == 20) return kPdgSigma0;       //  3212 / Sigma0
  if(geant_code == 21) return kPdgSigmaM;       //  3112 / Sigma-
  if(geant_code == 29) return kPdgAntiSigmaP;   // -3112 / \bar{Sigma+}
  if(geant_code == 28) return kPdgAntiSigma0;   // -3212 / \bar{Sigma0}
  if(geant_code == 27) return kPdgAntiSigmaM;   // -3112 / \bar{Sigma-}
  if(geant_code == 22) return kPdgXi0;          //  3322 / Xi0
  if(geant_code == 23) return kPdgXiM;          //  3312 / Xi-
  if(geant_code == 30) return kPdgAntiXi0;      // -3322 / \bar{Xi0}
  if(geant_code == 31) return kPdgAntiXiP;      // -3312 / \bar{Xi+}
  if(geant_code == 24) return kPdgOmegaM;       //  3334 / Omega-
  if(geant_code == 32) return kPdgAntiOmegaP;   // -3334 / \bar{Omega+}

  // some rare Geant3 codes that don't really need definitions in PDGCodes.h
  const int kPdgDeuteron = 1000010020; // pdg::IonPdgCode(2,1);
  const int kPdgTritium  = 1000010030; // pdg::IonPdgCode(3,1);
  const int kPdgAlpha    = 1000020040; // pdg::IonPdgCode(4,2);
  const int kPdgHe3      = 1000020030; // pdg::IonPdgCode(3,2);
  if(geant_code == 45) return kPdgDeuteron;
  if(geant_code == 46) return kPdgTritium;
  if(geant_code == 47) return kPdgAlpha;
  if(geant_code == 49) return kPdgHe3;

  LOG("PDG", pWARN)
        << "Can not convert geant code: " << geant_code << " to PDG";
  return 0;
}
//____________________________________________________________________________

