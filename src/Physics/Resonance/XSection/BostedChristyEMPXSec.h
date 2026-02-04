//____________________________________________________________________________
/*!

\class    genie::BostedChristyEMPXSec

\brief     Fit to inelastic cross sections for A(e,e')X
           valid for all W<3 GeV and all Q2<10 GeV2
          
\author   Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research
          based on fortran code provided on Peter Bosted's site: 
          https://userweb.jlab.org/~bosted/fits.html

\ref      1. M.E. Christy, P.E.Bosted, "Empirical fit to precision inclusive 
          electron-proton cross sections in the resonance region", PRC 81 (2010) 055213
          2. P.E.Bosted, M.E.Christy,  "Empirical fit to inelastic electron-deuteron 
          and electron-neutron resonance region transverse cross", PRC 77 (2008) 065206
          3. C. Maieron, T. W. Donnelly, and I. Sick, "Extended superscaling of electron 
          scattering from nuclei", PRC 65 (2001) 025502 


\created  April 3, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _BOOSTED_CHRISTY_EM_PXSEC_H_
#define _BOOSTED_CHRISTY_EM_PXSEC_H_

#include <array>
#include <map>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/NuclearState/FermiMomentumTable.h"

class XSecIntegratorI;

namespace genie {

class BostedChristyEMPXSec : public XSecAlgorithmI {

public:
  BostedChristyEMPXSec();
  BostedChristyEMPXSec(string config);
  virtual ~BostedChristyEMPXSec();

  // implement the XSecAlgorithmI interface
  double XSec           (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral       (const Interaction * i) const;
  bool   ValidProcess   (const Interaction * i) const;
  bool   ValidKinematics(const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);
  double sigmaR(int, double, double, bool) const;
  double sigmaNR(int, double, double, bool) const;
  void BranchingRatios(int, double&, double&) const;
  void FermiSmearingD(double, double, double&, double&, double&, double&, bool) const;
  void FermiSmearingA(double, double, double, double, double&, double&, double&, double&) const;
  double FitEMC(double, int) const;
  double MEC2009(int, double, double) const;

  bool   fUseMEC;                                      ///< account for MEC contribution?
  double fPM;                                          ///< mass parameter
  double fMP;                                          ///< mass parameter
  double fAM;                                          ///< mass parameter
  double fMD;                                          ///< deuterium mass
  double fMpi0;                                        ///< pion mass
  double fMeta;                                        ///< eta mass
  double fWmin;                                        ///< minimal W
  double fWmax;                                        ///< maximal W
  double fQ2min;                                       ///< minimal Q2
  double fQ2max;                                       ///< maximal Q2
  
  std::array<std::array<double, 3>, 7> fBRp;           ///< branching ratios of resonances for proton fit
  std::array<std::array<double, 3>, 7> fBRD;           ///< branching ratios of resonances for deterium fit
  
  std::array<int, 7> fAngRes;                          ///< resonance angular momentum
  
  std::array<double, 7> fMassRes;                      ///< resonance mass
  
  std::array<double, 7> fWidthRes;                     ///< resonance width
  
  std::array<std::array<double, 4>, 7> fRescoefTp;     ///< tunable parameters from Ref.1, Table III for resonance \sigma_T
  std::array<std::array<double, 4>, 7> fRescoefTD;     ///< tunable parameters from Ref.2, Table III for resonance \sigma_T
  std::array<std::array<double, 3>, 7> fRescoefL;      ///< tunable parameters from Ref.1, Table III for resonance \sigma_L
  
  std::array<std::array<double, 5>, 2> fNRcoefTp;      ///< tunable parameters from Ref.1, Table III for nonres bkg \sigma_T
  std::array<std::array<double, 5>, 2> fNRcoefTD;      ///< tunable parameters from Ref.1, Table IV  for nonres bkg \sigma_T
  std::array<double, 6>  fNRcoefL;                     ///< tunable parameters from Ref.1, Table III for nonres bkg \sigma_L
  std::array<double, 6>  fMECcoef;                     ///< tunable parameters for Eqs.(20), (21) Ref.2
  std::array<double, 8>  fMEC2009coef;                 ///< tunable parameters for MEC2009 function
  std::array<double, 13> fAfitcoef;                    ///< tunable parameters for nuclei fit
  
  std::array<double, 9> fEMCalpha;                    ///< tunable parameters for EMC fit
  std::array<double, 3> fEMCc;                        ///< tunable parameters for EMC fit
  
  map<int, double> fMEC2009p18;                            
  map<int, double> fKFTable; 
  map<int, double> fNucRmvE;

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace




#endif  // _BOOSTED_CHRISTY_EM_PXSEC_H_
