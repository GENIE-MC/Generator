//____________________________________________________________________________
/*!

\namespace genie::pdg

\brief     Utilities for improving the code readability when using PDG codes.
           E.g. a ..............................if( pdg::IsProton(pdg_code) )
           is much easier to read/check than....if( pdgc_code == 2212 )

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 06, 2004

\cpright   Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
           All rights reserved.
           For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PDG_UTILS_H_
#define _PDG_UTILS_H_

namespace genie {

namespace pdg
{
  int  IonPdgCodeToZ   (int ion_pdgc);
  int  IonPdgCodeToA   (int ion_pdgc);
  int  IonPdgCode      (int A, int Z);
  bool IsIon           (int pdgc);
  
  bool IsNeutrino      (int pdgc);
  bool IsAntiNeutrino  (int pdgc);
  bool IsNeutralLepton (int pdgc);
  bool IsChargedLepton (int pdgc);

  bool IsNuE           (int pdgc);
  bool IsNuMu          (int pdgc);
  bool IsNuTau         (int pdgc);
  bool IsAntiNuE       (int pdgc);
  bool IsAntiNuMu      (int pdgc);
  bool IsAntiNuTau     (int pdgc);
  
  bool IsElectron      (int pdgc);
  bool IsPositron      (int pdgc);
  bool IsMuon          (int pdgc);
  bool IsAntiMuon      (int pdgc);
  bool IsTau           (int pdgc);
  bool IsAntiTau       (int pdgc);
  
  bool IsQuark         (int pdgc);
  bool IsUQuark        (int pdgc);
  bool IsDQuark        (int pdgc);
  bool IsSQuark        (int pdgc);
  bool IsCQuark        (int pdgc);
  bool IsAntiQuark     (int pdgc);
  bool IsAntiUQuark    (int pdgc);
  bool IsAntiDQuark    (int pdgc);
  bool IsAntiSQuark    (int pdgc);
  bool IsAntiCQuark    (int pdgc);
  
  bool IsProton        (int pdgc);
  bool IsNeutron       (int pdgc);
  
  bool IsNeutronOrProton      (int pdgc);
  int  SwitchProtonNeutron    (int pdgc);
  int  Neutrino2ChargedLepton (int pdgc);

  bool IsBaryonResonance      (int pdgc);

}      // pdg namespace

}      // genie namespace

#endif // _PDG_UTILS_H_
