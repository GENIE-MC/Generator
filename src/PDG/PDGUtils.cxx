//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 02, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include "BaryonResonance/BaryonResUtils.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"

using namespace genie;

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
bool genie::pdg::IsIon(int pdgc)
{
  return (pdgc > 1000000000 && pdgc < 1999999999);
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
bool genie::pdg::IsNeutralLepton(int pdgc)
{
  bool is_neutral_lepton = IsNeutrino(pdgc) || IsAntiNeutrino(pdgc);
  return is_neutral_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsChargedLepton(int pdgc)
{
  bool is_neg_lepton = (pdgc ==  kPdgElectron) ||
                       (pdgc ==  kPdgMuon)     ||
                       (pdgc ==  kPdgTau);

  bool is_pos_lepton = (pdgc ==  kPdgPositron) ||
                       (pdgc ==  kPdgAntiMuon) ||
                       (pdgc ==  kPdgAntiTau);

  bool is_charged_lepton = is_neg_lepton || is_pos_lepton;
  return is_charged_lepton;
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

